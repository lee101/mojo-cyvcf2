"""A fast, compatible subset of :mod:`cyvcf2` for plain-text VCF files."""
from .vcf import VCF, Variant

Reader = VCF
VCFReader = VCF
__all__ = ["VCF", "Reader", "VCFReader", "Variant"]
__version__ = "0.1.0"
