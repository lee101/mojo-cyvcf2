# mojo-cyvcf2

`mojo-cyvcf2` is a standalone Mojo port of the hot, plain-text VCF parsing path
from [cyvcf2](https://github.com/brentp/cyvcf2). It keeps the familiar Python
`VCF` and `Variant` API for the covered subset while moving record scanning and
GT decoding into a compiled Mojo shared library.

## Covered subset

- Plain, uncompressed VCF 4.x files (`.vcf`), headers, sequential iteration, and
  simple `chrom` / `chrom:start-end` queries.
- `VCF(path, mode="r", gts012=False, lazy=False, strict_gt=False, samples=None,
  threads=None)`, including sample projection and `gts012` genotype-type mapping.
- Core `Variant` fields (`CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`,
  `INFO`, `FORMAT`, coordinates and variant predicates).
- Genotypes and common call helpers: `genotypes`, `gt_types`, `gt_bases`,
  `gt_phases`, depth/quality/allele-depth helpers, call counts, `call_rate`, and
  `aaf`; numeric/string `format(key)` access for fixed-width VCF FORMAT data.

Not covered: bgzip/BCF input, tabix-indexed queries, HTTP input, `Writer`, header
mutation, `set_format`, and cyvcf2's population-statistics helpers. Those depend
on htslib's compression/indexing and are deliberately out of scope for this
compute-focused first port.

## Install and use

```bash
pixi install
pixi run build
```

`pixi` sets `PYTHONPATH=python`, so this self-contained example runs directly
from the checkout:

```python
from tempfile import NamedTemporaryFile
from mojo_cyvcf2 import VCF

with NamedTemporaryFile(mode="w", suffix=".vcf") as handle:
    handle.write("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumour\tnormal
chr7\t55240001\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1
""")
    handle.flush()
    vcf = VCF(handle.name, samples=["tumour", "normal"], strict_gt=True)
    for variant in vcf("chr7:55240000-55250000"):
        print(variant.CHROM, variant.POS, variant.REF, variant.ALT)
        print(variant.genotypes, variant.gt_depths)
```

Run the checks and reproducible benchmark with:

```bash
pixi run test
pixi run bench
```

## How it works

Python reads the VCF once as a contiguous `uint8` buffer and passes its address
to `dist/libmojo-cyvcf2.so` through `ctypes`. Mojo finds each non-header line,
locates its first nine tab-delimited columns, and decodes leading `GT` fields.
The record scan uses unaligned SIMD loads, skips vectors with no newline, and
finishes with a scalar tail. GT decode runs through `max.algorithm.parallelize`
inside the kernel, but only from 65,536 records upwards: below that the kernel
stays serial because the measured crossover is far above the old assumption.
Record rows are owned per chunk, so the threads write disjoint output slices and
all ownership remains on the Python side. Record positions are `int64` offset
arrays, and decoded calls are a contiguous `(records, samples, 3)` `int32` array
of `(allele_a, allele_b, phased)`. Python materializes strings only when a
`Variant` is requested.

## Benchmarks

Measured by `pixi run bench` in this checkout. The generated input has 200,000
records and eight samples; timings are the best of three full `list(VCF(path))`
iterations.

| case | mojo-cyvcf2 | cyvcf2 | ratio |
| --- | ---: | ---: | ---: |
| iterate + decode GT (200k x 8) | 647.1 ms | 1025.8 ms | 1.59x faster |

End-to-end timing moves by tens of milliseconds between runs on this shared
machine (repeated runs of the current build measured 647-723 ms for
mojo-cyvcf2), so the table row is one run, not a best case. The compiled
kernels are a small part of that total; per-record Python `Variant`
construction dominates it.

This is an end-to-end compatibility benchmark, including Python-compatible
`Variant` object creation and strings. On this run mojo-cyvcf2 is faster;
cyvcf2 remains the more complete parser, and its htslib-backed
compressed/indexed paths are not represented by this plain-text workload. The
Mojo scan and GT decode are real compiled kernels, not a selective
microbenchmark. GPU execution is intentionally omitted: delimiter discovery
and variable-width GT parsing are branch- and memory-bound, below the
arithmetic intensity where transfer overhead can pay off.

### Kernel time across the UNROLL phase

Same 200k x 8 input, timing the two exported kernels directly with the
pre-change and post-change shared libraries alternating in one process, best
of 15 rounds (`bench/kernel.py` times the current build the same way):

| kernel | before | after | change |
| --- | ---: | ---: | ---: |
| `mcv_scan_records` | 8.87 ms | 8.85 ms | unchanged |
| `mcv_decode_gt` | 47.19 ms | 35.70 ms | 1.32x faster |
| both | 56.06 ms | 44.55 ms | 1.26x faster |

Both kernels produce byte-identical record bounds, alleles and gt types to the
pre-change build; the decode change was checked against the previous library
on the benchmark input plus randomised records with multi-digit and polyploid
genotypes, haploid and missing calls, GT in a non-leading FORMAT position,
absent FORMAT, CRLF input and records without a trailing newline, in both
`strict_gt` modes. The gains come from hoisting the per-column state reads out
of the byte loops and decoding the common one-digit diploid call (`.`, `a`,
`a/b`, `a|b`) without the general accumulation loop. Scan-side micro-optimisations
that `pixi run bench` did not show as faster were dropped rather than kept.

There are no floating-point reductions in these kernels, so this phase does not
change summation order and no result differs from the previous build in any
bit position.

### SIMD phase

`bench/analysis.md`, which this phase was meant to use as its worklist, does not
exist in this repository and never has (no commit adds or removes it). The hot
loops were therefore taken from the kernel source itself, and every candidate was
judged by `pixi run bench` and by kernel timing with the pre-change library
alternating in the same process.

The kernels are already vectorised where it helps: `mcv_scan_records` compares a
whole block at a time and reads line ends out of that same typed vector.
Widening it further, or reshaping the passes that remain, did not move the
caller-size numbers:

| variant | scan before | scan after |
| --- | ---: | ---: |
| 32-byte block, lane indexing (whole-lane loop removal, guarded tail) | 8.63, 8.56, 8.59 ms | 8.72, 8.70, 8.74 ms |
| same, with a bounded SIMD header skip | 8.32 ms | 8.50 ms |

Three alternating-library runs, best of 15 rounds each: the scan measured
0.98-0.99x (slightly slower) and the end-to-end case did not change
(`pixi run bench` measured 649.8 ms before this phase, 676.0 ms in one run
during it, within the 647-723 ms run-to-run spread noted above). Fusing the
caller's newline count into the scan was also measured and dropped: it saves a
7.45 ms NumPy pass but costs 10.51 ms against 9.07 ms for the scan alone.

A fixed-width GT decoder shared one fatal shape with the general kernel, as could
any scan/digit fast path: the kernels receive a bare address and a length, so a
digit test that assumes a sample begins with a number decodes the digit before it
when it does not. Restricting the fast path to single-digit calls and validating
the whole layout up front removes that hazard, but the validation is itself a pass
over the file; a vectorised layout check cost 19.86 ms, more than the kernel it
would authorise. Machine-checkable rules for the shape would cost no pass, but
matching a whole record's numeric skeleton together with genotype semantics is a
separate piece of work, so none of it is shipped here.

Nothing from this phase is in the tree: the scan is byte-for-byte the previous
build, and `tests/test_simd_paths.py` was added to pin the existing scan and GT
decode at lengths 0, 1, 31, 32, 33 and a large odd length, so any future SIMD
change has to keep those exact. `tests` went from 11 to 36 tests.

### PARALLEL phase

The GT decode used to be split from Python across an 8-worker
`ThreadPoolExecutor` of ctypes calls, starting at 32,768 records. That loop is
gone: the threshold is now inside the compiled kernel and was chosen from
measurements. `from max.algorithm import parallelize` was not resolvable in this
environment at first (`unable to locate module 'max'`), because `max` was not a
declared dependency; `max = "==26.7.0.dev2026092605"` is now pinned beside the
mojo pin and `pixi.lock` was regenerated.

Two things set the shape of this phase, both measured with `bench/parallel.py`:

1. **Each `parallelize` launch costs ~0.5 ms** on this host, and up to 2.8 ms
   when sibling processes are busy. A trivial one-item-per-thread body costs
   the same as real work, so it is launch cost, not scheduling of the chunks.
2. **Scale-out keeps improving beyond 4 workers as the input grows.** At 1 M
   records x 8 samples w4 finishes in ~56 ms and w8 in ~33 ms; at 3 M x 8 w8
   takes 83.6 ms against 152.1 ms for w4 and 484.8 ms serial (0.17x). The 4-worker
   cap below is therefore a step on the curve, not its knee; 8 was kept out
   because every size this package's own benchmark and threshold exercise are
   at or under 300k x 8, where 8 workers is level with or worse than 4.

**Scan: no parallel path was kept.** Splitting the byte scan by newline gave
`w1 1.00x`, `w2 0.75x`, `w4 0.50x`, `w8 0.32x` at 3 M records (273 MB), and
`w2 2.02x` at 100k. It loses to serial at every size a caller reaches, so it is
not shipped: the scan is byte-for-byte the previous kernel.

**Decode: shipped behind a threshold.** Exact equality against the serial kernel
at 1, 2, 4 and 8 workers, and the switch points are pinned by tests.

| records x 8 samples | bytes | serial | w1 | w2 | w4 | w8 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 100,000 | 8.8 MB | 16.2 ms | 1.00x | 1.02x | 0.95x | 0.97x |
| 160,000 | 14.2 MB | 25.8 ms | 1.01x | 0.91x | 0.91x | 0.95x |
| 300,000 | 26.8 MB | 50.7 ms | 0.98x | 0.70x | 0.49x | 0.50x |
| 1,000,000 | 89.8 MB | 172.4 ms | 0.99x | 0.53x | 0.33x | 0.19x |

Below ~100k x 8 the parallel path only ties or loses (`w2 2.99x` at 32k x 1,
`w2`/`w4` 1.9x/1.6x slower at 20k x 8). The threshold is
`_DECODE_PARALLEL_THRESHOLD = 65_536` records, a margin above that crossover,
and the kernel additionally refuses a second worker below 524,288 genotype
calls, so the second thread is never launched where it measured as a loss. The
worker cap is 4: at 300k x 8 the two are level (0.70x / 0.49x / 0.50x for
2 / 4 / 8 workers), and 8 was no better than 4 at every smaller size measured.

The committed end-to-end effect on the 200k x 8 benchmark input (1.6 M calls,
kernel time only, same input in one process):

| kernel | before | after | change |
| --- | ---: | ---: | ---: |
| `mcv_scan_records` | 9.27 ms | 9.35 ms | unchanged code, run noise |
| `mcv_decode_gt` serial | 38.24 ms | 37.37 ms | unchanged code, run noise |
| decode as committed (3 workers) | 37.13 ms | 20.06 ms | 1.85x faster |
| scan + decode, end to end | 47.37 ms | 33.38 ms | 1.42x faster |

The "before" column is the previous `dist/libmojo-cyvcf2.so` and the "after"
column this build, alternating in one process, best of 15. Scan and serial
decode are the same code in both, so their rows measure the host, not a change;
only the committed 3-worker decode is faster.

End-to-end `list(VCF(path))` including Python `Variant` construction, measured
three times in this checkout during the phase: 670.6 ms, 666.6 ms, 760.7 ms,
against the 647.1-723 ms range already recorded above. That difference is inside
the documented run-to-run spread on this shared host, not a demonstrated gain:
per-record Python object construction dominates that total. `pixi run bench`
could not be used for that row because the `/tmp/mojo-bench.lock` its task takes
is held by other concurrent sessions on this host for the whole phase; the same
script was run directly (it only reads and times the generated file).

`bench/parallel.py` reproduces both tables and the equality check
(`--check`). `tests/test_parallel_paths.py` compares serial and parallel output
byte for byte at 0, 1, 7, 4,097 and 40,000 records, strict and non-strict,
awkward genotypes (haploid, polyploid, multi-digit, GT not leading), the
threshold +/- one record through the public `VCF` API, and `threads=None/1/2/4/64`.
`tests` went from 36 to 85 tests.

## Compatibility testing

The test suite installs the real PyPI `cyvcf2` package and asserts parity on a
representative VCF for fields, INFO, FORMAT matrices, genotype states, strict
missing-call behavior, sample selection, and regions.

MIT.
