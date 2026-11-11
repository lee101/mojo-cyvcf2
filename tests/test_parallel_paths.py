"""Serial and parallel decode must produce identical output on the same input.

The parallel path is gated by `min_records` and by calls-per-worker inside the
kernel, so these tests pin both the equality and the switch points. Every check
compares serial output byte for byte, not within a tolerance.
"""
import numpy as np
import pytest

from mojo_cyvcf2._lib import addr, lib

SAMPLES = 8


def _vcf(records, samples=SAMPLES, calls=("0/0:12", "0|1:20", "1/1:9", "./.:.")):
    header = (
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(f"S{i}" for i in range(samples)) + "\n"
    )
    row_calls = "\t".join(calls[i % len(calls)] for i in range(samples))
    body = "".join(f"1\t{i + 1}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t{row_calls}\n" for i in range(records))
    return (header + body).encode()


def _scan(src):
    capacity = int(np.count_nonzero(src == 10)) + 1
    starts = np.empty(capacity, dtype=np.int64)
    ends = np.empty(capacity, dtype=np.int64)
    count = lib().mcv_scan_records(addr(src), len(src), addr(starts), addr(ends), capacity)
    return count, starts[:count], ends[:count]


def _decode(src, starts, ends, count, samples, strict, workers=None, min_records=0):
    total = count * samples
    alleles = np.full(max(1, total * 3), -7, dtype=np.int32)
    types = np.full(max(1, total), -7, dtype=np.int32)
    if workers is None:
        lib().mcv_decode_gt(addr(src), addr(starts), addr(ends), count, samples, strict,
                            addr(alleles), addr(types))
    else:
        lib().mcv_decode_gt_parallel(addr(src), addr(starts), addr(ends), count, samples, strict,
                                     addr(alleles), addr(types), workers, min_records)
    return alleles[:total * 3], types[:total]


@pytest.mark.parametrize("workers", [1, 2, 3, 4])
@pytest.mark.parametrize("records", [0, 1, 7, 4_097, 40_000])
@pytest.mark.parametrize("strict", [0, 1])
def test_parallel_decode_matches_serial(records, workers, strict):
    src = np.frombuffer(_vcf(records), dtype=np.uint8)
    count, starts, ends = _scan(src)
    assert count == records
    ref_alleles, ref_types = _decode(src, starts, ends, count, SAMPLES, strict)
    alleles, types = _decode(src, starts, ends, count, SAMPLES, strict, workers)
    np.testing.assert_array_equal(alleles, ref_alleles)
    np.testing.assert_array_equal(types, ref_types)


def test_parallel_decode_handles_awkward_genotypes():
    calls = (".", "0", "1", "./.", ".|1", "10/2", "1/2/3", ".:.", "0|0:5", "./.:.")
    for records in (5_000, 40_001):
        src = np.frombuffer(_vcf(records, samples=len(calls), calls=calls), dtype=np.uint8)
        count, starts, ends = _scan(src)
        ref_alleles, ref_types = _decode(src, starts, ends, count, len(calls), 0)
        for workers in (2, 4):
            alleles, types = _decode(src, starts, ends, count, len(calls), 0, workers)
            np.testing.assert_array_equal(alleles, ref_alleles)
            np.testing.assert_array_equal(types, ref_types)


def test_kernel_falls_back_to_serial_below_min_records():
    """A caller below the threshold must get the serial result, row for row."""
    records = 1_000
    src = np.frombuffer(_vcf(records), dtype=np.uint8)
    count, starts, ends = _scan(src)
    ref_alleles, ref_types = _decode(src, starts, ends, count, SAMPLES, 0)
    for workers in (2, 4):
        alleles, types = _decode(src, starts, ends, count, SAMPLES, 0, workers,
                                 min_records=records + 1)
        np.testing.assert_array_equal(alleles, ref_alleles)
        np.testing.assert_array_equal(types, ref_types)


def test_kernel_caps_workers_by_call_count():
    """Under 2 x 524,288 calls the kernel runs serial even when asked for 4."""
    records = 80_000
    src = np.frombuffer(_vcf(records), dtype=np.uint8)
    count, starts, ends = _scan(src)
    ref_alleles, ref_types = _decode(src, starts, ends, count, SAMPLES, 0)
    # 80,000 x 8 = 640,000 calls: under two CALLS_PER_WORKER units, so serial.
    alleles, types = _decode(src, starts, ends, count, SAMPLES, 0, 4, min_records=0)
    np.testing.assert_array_equal(alleles, ref_alleles)
    np.testing.assert_array_equal(types, ref_types)


def test_vcf_reader_parallel_and_serial_agree_at_threshold_boundaries(tmp_path):
    """The committed threshold, +/- one record, must give the serial kernel's output."""
    from mojo_cyvcf2.vcf import _DECODE_PARALLEL_THRESHOLD, VCF

    for records in (_DECODE_PARALLEL_THRESHOLD - 1, _DECODE_PARALLEL_THRESHOLD,
                    _DECODE_PARALLEL_THRESHOLD + 1, 120_000):
        path = tmp_path / f"n{records}.vcf"
        path.write_bytes(_vcf(records))
        vcf = VCF(path, threads=4)
        src = np.frombuffer(path.read_bytes(), dtype=np.uint8)
        count, starts, ends = _scan(src)
        assert count == records
        ref_alleles, ref_types = _decode(src, starts, ends, count, SAMPLES, 0)
        np.testing.assert_array_equal(vcf._alleles.reshape(-1), ref_alleles)
        np.testing.assert_array_equal(vcf._types.reshape(-1), ref_types)


@pytest.mark.parametrize("threads", [None, 1, 2, 4, 64])
def test_vcf_reader_matches_reference_above_threshold(tmp_path, threads):
    from mojo_cyvcf2 import VCF

    path = tmp_path / "many.vcf"
    path.write_bytes(_vcf(40_000))
    vcf = VCF(path, threads=threads)
    assert len(vcf._starts) == 40_000
    # Decode the same file through the plain serial kernel and compare.
    src = np.frombuffer(path.read_bytes(), dtype=np.uint8)
    count, starts, ends = _scan(src)
    ref_alleles, ref_types = _decode(src, starts, ends, count, SAMPLES, 0)
    np.testing.assert_array_equal(vcf._alleles.reshape(-1), ref_alleles)
    np.testing.assert_array_equal(vcf._types.reshape(-1), ref_types)
    # `threads` caps the workers at 4; None means the default cap.
    assert vcf.threads == (4 if threads is None else min(threads, 4))
