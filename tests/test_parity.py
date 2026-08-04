"""Behavioural parity checks against cyvcf2 on published VCF conventions."""
from __future__ import annotations

import numpy as np
import pytest

import mojo_cyvcf2 as mojo

cyvcf2 = pytest.importorskip("cyvcf2")

VCF_TEXT = """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description=Allele frequency>
##INFO=<ID=DB,Number=0,Type=Flag,Description=Known>
##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=Allele depths>
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=Quality>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
1\t101\trs1\tA\tG\t50\tPASS\tAF=0.25;DB\tGT:DP:AD:GQ\t0/0:12:12,0:99\t0|1:20:10,10:48\t1/1:9:0,9:12
1\t205\t.\tAT\tA\t.\tq10\t.\tGT:DP:AD\t0/1:8:4,4\t./.:.:.\t0/0:7:7,0
2\t9\trs3\tC\tT,G\t19.5\t.\tAF=0.2,0.1\tGT:DP\t1/2:15\t0/2:11\t0/.:3
"""


@pytest.fixture()
def vcf_path(tmp_path):
    path = tmp_path / "calls.vcf"
    path.write_text(VCF_TEXT)
    return path


def test_core_fields_and_iteration_match_cyvcf2(vcf_path):
    ours, upstream = list(mojo.VCF(vcf_path)), list(cyvcf2.VCF(str(vcf_path)))
    assert len(ours) == len(upstream) == 3
    for a, b in zip(ours, upstream):
        assert (a.CHROM, a.POS, a.ID, a.REF, a.ALT) == (b.CHROM, b.POS, b.ID, b.REF, b.ALT)
        assert a.QUAL == b.QUAL
        assert a.FILTER == b.FILTER
        assert a.start == b.start and a.end == b.end
        assert a.is_snp == b.is_snp and a.is_indel == b.is_indel and a.is_deletion == b.is_deletion


def test_info_format_and_genotypes_match_cyvcf2(vcf_path):
    ours, upstream = list(mojo.VCF(vcf_path)), list(cyvcf2.VCF(str(vcf_path)))
    for a, b in zip(ours, upstream):
        assert a.genotypes == b.genotypes
        assert np.array_equal(a.gt_types, b.gt_types)
        assert np.array_equal(a.gt_phases, b.gt_phases)
        assert np.array_equal(a.gt_depths, b.gt_depths)
        assert np.array_equal(a.gt_ref_depths, b.gt_ref_depths)
        assert np.array_equal(a.gt_alt_depths, b.gt_alt_depths)
        assert np.array_equal(a.gt_quals, b.gt_quals)
        for key in a.FORMAT:
            if key != "GT":
                assert np.array_equal(a.format(key), b.format(key), equal_nan=True)
    assert ours[0].INFO["AF"] == pytest.approx(upstream[0].INFO["AF"])
    assert ours[0].INFO["DB"] is True
    assert ours[2].INFO["AF"] == pytest.approx(tuple(upstream[2].INFO["AF"]))


def test_gts012_strict_gt_sample_selection_and_region(vcf_path):
    ours = list(mojo.VCF(vcf_path, gts012=True, strict_gt=True, samples=["S2", "S3"]))
    upstream = list(cyvcf2.VCF(str(vcf_path), gts012=True, strict_gt=True, samples=["S2", "S3"]))
    for a, b in zip(ours, upstream):
        assert a.genotypes == b.genotypes
        assert np.array_equal(a.gt_types, b.gt_types)
    reader = mojo.VCF(vcf_path)
    assert [(v.CHROM, v.POS) for v in reader("1:150-250")] == [("1", 205)]
    assert [(v.CHROM, v.POS) for v in reader("2")] == [("2", 9)]
    selected = ["S1", "S2", "S3"]
    reader = mojo.VCF(vcf_path, samples=selected)
    ours = next(iter(reader))
    upstream = next(iter(cyvcf2.VCF(str(vcf_path), samples=selected)))
    assert reader.samples == selected
    assert ours.genotypes == upstream.genotypes


def test_header_and_unsupported_inputs(vcf_path, tmp_path):
    reader = mojo.VCF(vcf_path)
    assert reader.samples == ["S1", "S2", "S3"]
    assert "##fileformat=VCFv4.2" in reader.raw_header
    with pytest.raises(ValueError, match="plain-text"):
        mojo.VCF(tmp_path / "calls.vcf.gz")


def test_documented_constructor_compatibility_arguments(vcf_path):
    reader = mojo.VCF(vcf_path, mode="rb", lazy=True, threads=2)
    assert reader.lazy is True and reader.threads == 2
    with pytest.raises(ValueError, match="read mode"):
        mojo.VCF(vcf_path, mode="w")


def test_nonleading_gt_and_partial_or_haploid_calls_match_cyvcf2(tmp_path):
    path = tmp_path / "edge-calls.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=Depth>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n"
        "1\t1\t.\tA\tG,T\t.\t.\t.\tDP:GT\t4:0/.\t5:1\t6:./2\n"
    )
    for strict in (False, True):
        ours, upstream = next(iter(mojo.VCF(path, strict_gt=strict))), next(iter(cyvcf2.VCF(str(path), strict_gt=strict)))
        assert ours.genotypes == upstream.genotypes
        assert np.array_equal(ours.gt_types, upstream.gt_types)


def test_simd_scan_tail_and_compact_decoding(vcf_path):
    vcf_path.write_text(VCF_TEXT.rstrip("\n"))
    ours, upstream = list(mojo.VCF(vcf_path)), list(cyvcf2.VCF(str(vcf_path)))
    assert len(ours) == len(upstream) == 3
    assert ours[-1].POS == upstream[-1].POS == 9
    assert ours[0].gt_types.dtype == np.int32
    assert np.array_equal(ours[-1].gt_types, upstream[-1].gt_types)


def test_documented_variant_helpers_and_format_types(vcf_path):
    first, second, third = list(mojo.VCF(vcf_path))
    assert first.FILTERS == []
    assert second.FILTERS == ["q10"]
    assert first.var_type == "snp" and second.var_type == "indel"
    assert first.var_subtype == "ts" and second.var_subtype == "del"
    assert first.is_transition and not first.is_mnp and not first.is_sv
    assert np.array_equal(first.format("GT"), first.gt_bases)
    assert first.num_called == 3 and first.num_unknown == 0
    assert first.num_het == 1 and first.num_hom_ref == 1 and first.num_hom_alt == 1
    assert first.call_rate == 1.0 and first.aaf == pytest.approx(0.5)
    assert third.num_called == 3 and third.call_rate == 1.0

    structural = vcf_path.with_name("predicates.vcf")
    structural.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tAT\tGC\t.\t.\t.\n"
        "1\t2\t.\tA\t<DEL>\t.\t.\t.\n"
    )
    mnp, sv = list(mojo.VCF(structural))
    assert mnp.is_mnp and mnp.var_type == "mnp"
    assert sv.is_sv and sv.var_type == "sv"


def test_no_sample_records_crlf_tail_and_format_range_protection(tmp_path):
    path = tmp_path / "no-samples.vcf"
    path.write_bytes(
        b"##fileformat=VCFv4.2\r\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\r\n"
        b"1\t7\t.\tA\tG\t.\t.\t.\tGT\r"
    )
    record = next(iter(mojo.VCF(path)))
    assert record.POS == 7 and record.FORMAT == ["GT"] and record.genotypes == []

    rich = tmp_path / "format-types.vcf"
    rich.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
        "1\t1\t.\tA\tG\t.\t.\t.\tGT:GL:XX:DP\t0/1:-1.5,0.0:abc:10\t0/0:-2.0,1.0:def:1\n"
    )
    record = next(iter(mojo.VCF(rich)))
    assert record.gt_depths.tolist() == [-1, -1]
    assert record.format("GL").dtype == np.float64
    assert record.format("XX").tolist() == [["abc"], ["def"]]
    overflow = rich.with_name("overflow.vcf")
    overflow.write_text(rich.read_text().replace("abc:10", "abc:2147483648"))
    with pytest.raises(OverflowError, match="int32"):
        next(iter(mojo.VCF(overflow))).format("DP")
