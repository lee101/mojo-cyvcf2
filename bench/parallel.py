"""Serial vs parallel GT-decode kernel: equality check and crossover timing.

Both entry points run the same per-record code, so the comparison is exact
(byte-identical arrays), never a tolerance. The timing column adds the
kernel's own gate on calls per worker, so the numbers here are raw worker counts
at or above any size the gate lets through.

    pixi run python bench/parallel.py --check   # equality at 1, 2, 4, 8 workers
    pixi run python bench/parallel.py           # timing table and crossover
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, ROOT)

from mojo_cyvcf2._lib import addr, lib  # noqa: E402

BENCH_SAMPLES = 8
WORKERS = (1, 2, 4, 8)
REPEAT = 9


def make_input(records: int, samples: int = BENCH_SAMPLES) -> bytes:
    header = (
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(f"S{i}" for i in range(samples)) + "\n"
    )
    calls = "\t".join(("0/0:12", "0|1:20", "1/1:9", "./.:.")[i % 4] for i in range(samples))
    body = "".join(
        f"1\t{i + 1}\trs{i + 1}\tA\tG\t50\tPASS\t.\tGT:DP\t{calls}\n" for i in range(records)
    )
    return (header + body).encode()


def scan(src):
    n = len(src)
    capacity = int(np.count_nonzero(src == 10)) + 1
    starts = np.empty(capacity, dtype=np.int64)
    ends = np.empty(capacity, dtype=np.int64)
    count = lib().mcv_scan_records(addr(src), n, addr(starts), addr(ends), capacity)
    return count, starts[:count], ends[:count]


def best(fn, repeat: int = REPEAT) -> float:
    result = float("inf")
    for _ in range(repeat):
        start = time.perf_counter()
        fn()
        result = min(result, time.perf_counter() - start)
    return result


def decode_call(fn, src, starts, ends, count, samples, alleles, types, workers=None,
                min_records=0):
    """Bind every argument now; a late-bound loop variable would time the wrong count."""
    if workers is None:
        return lambda: fn(addr(src), addr(starts), addr(ends), count, samples, 0,
                          addr(alleles), addr(types))
    return lambda: fn(addr(src), addr(starts), addr(ends), count, samples, 0,
                      addr(alleles), addr(types), workers, min_records)


def outputs(src, starts, ends, count, samples, workers=None):
    total = count * samples
    alleles = np.full(max(1, total * 3), -7, dtype=np.int32)
    types = np.full(max(1, total), -7, dtype=np.int32)
    if workers is None:
        lib().mcv_decode_gt(addr(src), addr(starts), addr(ends), count, samples, 0,
                            addr(alleles), addr(types))
    else:
        lib().mcv_decode_gt_parallel(addr(src), addr(starts), addr(ends), count, samples, 0,
                                     addr(alleles), addr(types), workers, 0)
    return alleles, types


def check(src, starts, ends, count, samples) -> None:
    ref_alleles, ref_types = outputs(src, starts, ends, count, samples)
    for workers in WORKERS:
        alleles, types = outputs(src, starts, ends, count, samples, workers)
        np.testing.assert_array_equal(alleles, ref_alleles)
        np.testing.assert_array_equal(types, ref_types)


def timing(records: int, samples: int) -> list[tuple[int, float, float]]:
    src = np.frombuffer(make_input(records, samples), dtype=np.uint8)
    count, starts, ends = scan(src)
    total = count * samples
    alleles = np.empty(max(1, total * 3), dtype=np.int32)
    types = np.empty(max(1, total), dtype=np.int32)
    serial = best(decode_call(lib().mcv_decode_gt, src, starts, ends, count, samples, alleles, types))
    row = []
    for workers in WORKERS:
        t = best(decode_call(lib().mcv_decode_gt_parallel, src, starts, ends, count, samples,
                             alleles, types, workers))
        row.append((workers, serial, t))
    return row


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="equality only")
    parser.add_argument("--sizes", default="5000,20000,50000,100000,200000")
    parser.add_argument("--samples", type=int, default=BENCH_SAMPLES)
    args = parser.parse_args()
    sizes = [int(s) for s in args.sizes.split(",") if s]
    for records in sizes:
        src = np.frombuffer(make_input(records, args.samples), dtype=np.uint8)
        count, starts, ends = scan(src)
        check(src, starts, ends, count, args.samples)
        print(f"{count} records x {args.samples} samples: serial and parallel decode "
              f"identical at workers={WORKERS}", flush=True)
        if args.check:
            continue
        row = timing(records, args.samples)
        cells = [f"serial {row[0][1] * 1e3:7.2f} ms"]
        for workers, serial, t in row:
            cells.append(f"w{workers} {t * 1e3:7.2f} ms ({t / serial:4.2f}x)")
        print("   ".join(cells), flush=True)


if __name__ == "__main__":
    main()
