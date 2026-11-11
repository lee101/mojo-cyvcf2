"""Edge lengths for the SIMD scan so its scalar tail is always exercised.

The scanner's block width is 32 bytes. Lengths around that boundary, plus a large
odd length, are where a vector body and its remainder handling disagree; these
tests pin the record bounds the scanner reports and the GT decode built on them.
"""
from __future__ import annotations

import numpy as np
import pytest

import mojo_cyvcf2
from mojo_cyvcf2._lib import addr, lib

W = 32
LENGTHS = [0, 1, W - 1, W, W + 1, 1001]

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
)
ROW = "1\t{}\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1|1\n"


def _scan(tmp_path, text: str):
    """Run the exported scanner directly and return (count, starts, ends)."""
    data = np.frombuffer(text.encode(), dtype=np.uint8).copy()
    capacity = max(2, len(data) + 2)
    starts = np.zeros(capacity, dtype=np.int64)
    ends = np.zeros(capacity, dtype=np.int64)
    buffer = np.zeros(max(1, len(data)), dtype=np.uint8)
    buffer[: len(data)] = data
    count = lib().mcv_scan_records(addr(buffer), len(data), addr(starts), addr(ends), capacity)
    return count, starts[:count], ends[:count], buffer


@pytest.mark.parametrize("length", LENGTHS)
def test_scanner_finds_every_record_at_tail_lengths(tmp_path, length):
    text = HEADER + "".join(ROW.format(i + 1) for i in range(length))
    count, starts, ends, buffer = _scan(tmp_path, text)
    assert count == length
    lines = text.encode()[len(HEADER.encode()):].split(b"\n")[:-1]
    assert len(lines) == length
    assert [bytes(buffer[s:e]) for s, e in zip(starts, ends)] == lines


@pytest.mark.parametrize("length", LENGTHS)
def test_final_line_without_newline_counts_once(tmp_path, length):
    text = HEADER + "".join(ROW.format(i + 1) for i in range(length))
    if length:
        text = text[:-1]  # drop the last newline: the final record runs to the buffer end
    count, starts, ends, buffer = _scan(tmp_path, text)
    assert count == length
    if length:
        row = ROW.format(length).encode()
        assert bytes(buffer[starts[-1]:ends[-1]]) == row[:-1]


@pytest.mark.parametrize("length", LENGTHS)
def test_crlf_lines_at_tail_lengths(tmp_path, length):
    text = HEADER.replace("\n", "\r\n") + "".join(ROW.format(i + 1).replace("\n", "\r\n") for i in range(length))
    count, starts, ends, buffer = _scan(tmp_path, text)
    assert count == length
    row = ROW.format(1).rstrip("\n").encode()
    if length:
        assert bytes(buffer[starts[0]:ends[0]]) == row


@pytest.mark.parametrize("length", LENGTHS)
def test_gt_decode_at_tail_lengths(tmp_path, length):
    path = tmp_path / f"n{length}.vcf"
    path.write_text(HEADER + "".join(ROW.format(i + 1) for i in range(length)))
    reader = mojo_cyvcf2.VCF(path)
    assert len(list(reader)) == length
    expected = [[[0, 1, 0], [1, 1, 1]] for _ in range(length)]
    assert reader._alleles.tolist() == expected
    assert reader._types.tolist() == [[1, 3] for _ in range(length)]


def test_mixed_layout_still_matches_expected_decode(tmp_path):
    """The decode path and its fallback agree on a file mixing call shapes."""
    path = tmp_path / "mixed.vcf"
    calls = ["0/0", "0|1", "1/1", "./.", "1", "0", ".", "2/1", "0|.", "10/12", "0/1:5"]
    header = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(f"S{i}" for i in range(len(calls))) + "\n"
    )
    path.write_text(header + "1\t1\t.\tA\tG,T\t.\t.\t.\tGT\t" + "\t".join(calls) + "\n")
    record = next(iter(mojo_cyvcf2.VCF(path)))
    assert record.genotypes == [
        [0, 0, False], [0, 1, True], [1, 1, False], [-1, -1, False], [1, True],
        [0, True], [-1, True], [2, 1, False], [0, -1, True], [10, 12, False], [0, 1, False],
    ]
