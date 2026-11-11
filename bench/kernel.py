"""Kernel-only timing for the scan and GT decode hot path.

The end-to-end `pixi run bench` case also includes Python `Variant` object
creation, which is noisy and hides changes inside the compiled kernels.  This
script times the two exported kernels directly on the same 200k x 8 input.
"""
from __future__ import annotations

import sys
import time

import numpy as np

from bench import make_vcf
from mojo_cyvcf2._lib import addr, lib

REPEAT = 7


def best(fn) -> float:
    result = float("inf")
    for _ in range(REPEAT):
        start = time.perf_counter()
        fn()
        result = min(result, time.perf_counter() - start)
    return result


def main() -> None:
    path = make_vcf()
    with open(path, "rb") as handle:
        data = handle.read()
    src = np.frombuffer(data, dtype=np.uint8)
    n = len(src)
    capacity = int(np.count_nonzero(src == 10)) + 1
    starts = np.empty(capacity, dtype=np.int64)
    ends = np.empty(capacity, dtype=np.int64)
    scan, decode = lib().mcv_scan_records, lib().mcv_decode_gt

    count = scan(addr(src), n, addr(starts), addr(ends), capacity)
    starts, ends = starts[:count], ends[:count]
    samples = 8
    alleles = np.empty(count * samples * 3, dtype=np.int32)
    types = np.empty(count * samples, dtype=np.int32)

    t_scan = best(lambda: scan(addr(src), n, addr(starts), addr(ends), capacity))
    t_decode = best(lambda: decode(addr(src), addr(starts), addr(ends), count, samples, 0, addr(alleles), addr(types)))
    print(f"records={count} bytes={n}")
    print(f"scan   {t_scan * 1e3:8.2f} ms  {n / t_scan / 1e6:6.2f} GB/s")
    print(f"decode {t_decode * 1e3:8.2f} ms  {n / t_decode / 1e6:6.2f} GB/s")
    print(f"total  {(t_scan + t_decode) * 1e3:8.2f} ms")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
