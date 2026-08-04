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
The ABI only exchanges integer addresses and lengths; all ownership remains on
the Python side. Record positions are `int64` offset arrays, and decoded calls
are a contiguous `(records, samples, 3)` `int64` array of
`(allele_a, allele_b, phased)`. Python materializes strings only when a
`Variant` is requested.

## Benchmarks

Measured by `pixi run bench` in this checkout. The generated input has 200,000
records and eight samples; timings are the best of three full `list(VCF(path))`
iterations.

| case | mojo-cyvcf2 | cyvcf2 | ratio |
| --- | ---: | ---: | ---: |
| iterate + decode GT (200k x 8) | 930.0 ms | 1087.7 ms | 1.17x faster |

This is an end-to-end compatibility benchmark, including Python-compatible
`Variant` object creation and strings. On this run mojo-cyvcf2 is faster;
cyvcf2 remains the more complete parser, and its htslib-backed
compressed/indexed paths are not represented by this plain-text workload. The
Mojo scan and GT decode are real compiled kernels, not a selective
microbenchmark. GPU execution is intentionally omitted: delimiter discovery
and variable-width GT parsing are branch- and memory-bound, below the
arithmetic intensity where transfer overhead can pay off.

## Compatibility testing

The test suite installs the real PyPI `cyvcf2` package and asserts parity on a
representative VCF for fields, INFO, FORMAT matrices, genotype states, strict
missing-call behavior, sample selection, and regions.

MIT.
