"""Benchmark plain-text VCF scanning through the required pixi task."""
from __future__ import annotations

import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))

import cyvcf2
import mojo_cyvcf2


def make_vcf(records: int = 200_000, samples: int = 8) -> str:
    header = (
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n"
        "##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
        + "\t".join(f"S{i}" for i in range(samples)) + "\n"
    )
    calls = "\t".join(("0/0:12:12,0", "0|1:20:10,10", "1/1:9:0,9", "./.:.:.")[i % 4] for i in range(samples))
    row = "1\t{}\trs{}\tA\tG\t50\tPASS\tAF=0.25\tGT:DP:AD\t{}\n"
    path = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False).name
    with open(path, "w") as f:
        f.write(header)
        for i in range(records): f.write(row.format(i + 1, i + 1, calls))
    return path


def best(fn, repeat: int = 3) -> float:
    result = float("inf")
    for _ in range(repeat):
        start = time.perf_counter(); fn(); result = min(result, time.perf_counter() - start)
    return result


def main() -> None:
    path = make_vcf()
    try:
        mojo_cyvcf2.VCF(path)  # build/warm outside timing
        cases = [("iterate + decode GT (200k x 8)", lambda: list(mojo_cyvcf2.VCF(path)), lambda: list(cyvcf2.VCF(path)))]
        print("| case | mojo-cyvcf2 | cyvcf2 | ratio |")
        print("| --- | ---: | ---: | ---: |")
        for name, ours, theirs in cases:
            a, b = best(ours), best(theirs)
            print(f"| {name} | {a * 1e3:.1f} ms | {b * 1e3:.1f} ms | {b / a:.2f}x {'faster' if a < b else 'slower'} |")
    finally:
        os.unlink(path)


if __name__ == "__main__": main()
