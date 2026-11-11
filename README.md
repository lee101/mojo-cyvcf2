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
finishes with a scalar tail. For at least 32,768 records, independent GT decode
chunks run in a bounded CPU thread pool; smaller inputs stay serial. The ABI
only exchanges integer addresses and lengths, so the threads operate directly
on disjoint slices of the same NumPy buffers and all ownership remains on the
Python side. Record positions are `int64` offset arrays, and decoded calls are
a contiguous `(records, samples, 3)` `int32` array of
`(allele_a, allele_b, phased)`. Python materializes strings only when a
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

## Compatibility testing

The test suite installs the real PyPI `cyvcf2` package and asserts parity on a
representative VCF for fields, INFO, FORMAT matrices, genotype states, strict
missing-call behavior, sample selection, and regions.

MIT.
