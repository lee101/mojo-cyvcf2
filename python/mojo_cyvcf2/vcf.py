"""Compatibility objects built on the compiled byte scanner.

The reader deliberately targets uncompressed, tab-delimited VCF 4.x files.  The
Mojo layer owns record discovery and genotype decoding; Python only materializes
the strings requested by the public API.
"""
from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ._lib import addr, lib

HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


class INFO(dict):
    """cyvcf2-style mapping for a record's INFO column."""

    def get(self, key: str, default: Any = None) -> Any:  # type: ignore[override]
        return super().get(key, default)


def _coerce(value: str) -> Any:
    if value == ".":
        return None
    if "," in value:
        return tuple(_coerce(part) for part in value.split(","))
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def _parse_info(text: str) -> INFO:
    result = INFO()
    if text in ("", "."):
        return result
    for item in text.split(";"):
        key, sep, value = item.partition("=")
        result[key] = _coerce(value) if sep else True
    return result


def _int32_array(values: list[Any]) -> np.ndarray:
    """Make an ABI-compatible integer FORMAT array without truncation."""
    low, high = np.iinfo(np.int32).min, np.iinfo(np.int32).max
    if any(not isinstance(value, (int, np.integer)) or value < low or value > high for value in values):
        raise OverflowError("FORMAT integer value cannot be represented as int32")
    return np.asarray(values, dtype=np.int32)


def _decode_gt_text(text: str | None, strict: bool) -> tuple[int, int, int, int]:
    if text is None:
        return -1, -1, 0, UNKNOWN
    phased = int("|" in text)
    parts = text.replace("|", "/").split("/")
    values = []
    for part in parts[:2]:
        try: values.append(int(part))
        except ValueError: values.append(-1)
    if len(values) == 1: values.append(-1)
    first, second = values
    haploid = len(parts) == 1
    missing = first < 0 or (second < 0 and not haploid)
    if strict and missing: state = UNKNOWN
    elif missing and first < 0 and second < 0: state = UNKNOWN
    elif missing: state = HOM_REF if first == 0 or second == 0 else HET
    elif haploid: state = HOM_REF if first == 0 else HOM_ALT
    elif first == second: state = HOM_REF if first == 0 else HOM_ALT
    else: state = HET
    return first, second, phased, state


@dataclass(slots=True)
class Variant:
    """One VCF record, with cyvcf2's commonly used attributes and helpers."""

    _reader: "VCF"
    _index: int
    _fields: tuple[str, ...]

    @property
    def CHROM(self) -> str: return self._fields[0]

    @property
    def POS(self) -> int: return int(self._fields[1])

    @property
    def ID(self) -> str | None: return None if self._fields[2] == "." else self._fields[2]

    @property
    def REF(self) -> str: return self._fields[3]

    @property
    def ALT(self) -> list[str]: return self._fields[4].split(",") if self._fields[4] != "." else []

    @property
    def QUAL(self) -> float | None: return None if self._fields[5] == "." else float(self._fields[5])

    @property
    def FILTER(self) -> str | None: return None if self._fields[6] in (".", "PASS") else self._fields[6]

    @property
    def FILTERS(self) -> list[str]: return [] if self.FILTER is None else self.FILTER.split(";")

    @property
    def INFO(self) -> INFO: return _parse_info(self._fields[7])

    @property
    def FORMAT(self) -> list[str]: return [] if len(self._fields) < 9 else self._fields[8].split(":")

    @property
    def start(self) -> int: return self.POS - 1

    @property
    def end(self) -> int: return self.start + len(self.REF)

    @property
    def var_type(self) -> str:
        if self.is_snp: return "snp"
        if self.is_mnp: return "mnp"
        if self.is_sv: return "sv"
        if self.is_indel: return "indel"
        return "unknown"

    @property
    def var_subtype(self) -> str:
        if self.is_transition: return "ts"
        if self.is_snp: return "tv"
        if self.is_deletion: return "del"
        if self.is_indel: return "ins"
        return "unknown"

    @property
    def is_snp(self) -> bool:
        return len(self.REF) == 1 and bool(self.ALT) and all(len(a) == 1 and not a.startswith("<") for a in self.ALT)

    @property
    def is_mnp(self) -> bool:
        return len(self.REF) > 1 and bool(self.ALT) and all(len(a) == len(self.REF) for a in self.ALT)

    @property
    def is_indel(self) -> bool:
        return bool(self.ALT) and any(len(a) != len(self.REF) for a in self.ALT)

    @property
    def is_deletion(self) -> bool: return bool(self.ALT) and any(len(a) < len(self.REF) for a in self.ALT)

    @property
    def is_sv(self) -> bool: return bool(self.ALT) and any(a.startswith("<") or "]" in a or "[" in a for a in self.ALT)

    @property
    def is_transition(self) -> bool:
        return self.is_snp and all((self.REF, alt) in (("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")) for alt in self.ALT)

    @property
    def genotypes(self) -> list[list[int | bool]]:
        values = self._reader._alleles[self._index]
        gt_column = self.FORMAT.index("GT") if "GT" in self.FORMAT else -1
        result = []
        for (a, b, phased), sample in zip(values, self._fields[9:]):
            call = sample.split(":")[gt_column] if gt_column >= 0 and len(sample.split(":")) > gt_column else None
            if call is not None and "/" not in call and "|" not in call:
                result.append([int(a), True])
            else:
                result.append([int(a), int(b), bool(phased)])
        return result

    @property
    def gt_types(self) -> np.ndarray:
        values = self._reader._types[self._index].copy()
        if self._reader.gts012:
            unknown = values == UNKNOWN
            hom_alt = values == HOM_ALT
            values[unknown] = HOM_ALT
            values[hom_alt] = UNKNOWN
        return values

    @property
    def gt_phases(self) -> np.ndarray: return self._reader._alleles[self._index, :, 2].astype(bool)

    @property
    def gt_bases(self) -> np.ndarray:
        bases = []
        for a, b, phased in self.genotypes:
            separator = "|" if phased else "/"
            left = "." if a < 0 else self.REF if a == 0 else self.ALT[a - 1]
            right = "." if b < 0 else self.REF if b == 0 else self.ALT[b - 1]
            bases.append(f"{left}{separator}{right}")
        return np.asarray(bases, dtype=object)

    def _format_strings(self, key: str) -> list[str | None]:
        if key not in self.FORMAT:
            raise KeyError(key)
        column = self.FORMAT.index(key)
        return [sample.split(":")[column] if len(sample.split(":")) > column else None for sample in self._fields[9:]]

    def format(self, key: str) -> np.ndarray:
        """Return a FORMAT field as a NumPy vector/matrix, matching cyvcf2 use cases."""
        if key == "GT":
            return self.gt_bases
        vals = self._format_strings(key)
        if any(v is not None and "," in v for v in vals):
            width = max((len(v.split(",")) for v in vals if v not in (None, ".")), default=1)
            missing = -2147483648
            rows = [([missing + j for j in range(width)] if v in (None, ".") else [_coerce(x) for x in v.split(",")] + [missing + j for j in range(len(v.split(",")), width)]) for v in vals]
            flat = [value for row in rows for value in row]
            if all(isinstance(value, (int, np.integer)) for value in flat):
                return _int32_array(flat).reshape(len(rows), width)
            if all(isinstance(value, (int, float, np.number)) for value in flat):
                return np.asarray(rows, dtype=np.float64)
            return np.asarray(rows, dtype=object)
        if any(v is None or v == "." for v in vals):
            present = [_coerce(v) for v in vals if v not in (None, ".")]
            if all(isinstance(x, int) for x in present):
                return _int32_array([-2147483648 if v in (None, ".") else _coerce(v) for v in vals]).reshape(-1, 1)
            return np.asarray([np.nan if v in (None, ".") else _coerce(v) for v in vals], dtype=float).reshape(-1, 1)
        parsed = [_coerce(v) for v in vals]
        if all(isinstance(value, int) for value in parsed):
            return _int32_array(parsed).reshape(-1, 1)
        if all(isinstance(value, (int, float)) for value in parsed):
            return np.asarray(parsed, dtype=np.float64).reshape(-1, 1)
        return np.asarray(parsed, dtype=object).reshape(-1, 1)

    @property
    def gt_depths(self) -> np.ndarray:
        if "AD" not in self.FORMAT:
            return np.full(len(self._reader.samples), -1, dtype=np.int32)
        return self._depths("DP")
    @property
    def gt_quals(self) -> np.ndarray: return self._depths("GQ")

    def _depths(self, key: str) -> np.ndarray:
        try:
            values = self.format(key).astype(np.int32, copy=False).reshape(-1)
            values[values <= -2147483640] = -1
            return values
        except KeyError: return np.full(len(self._reader.samples), -1, dtype=np.int32)

    @property
    def gt_ref_depths(self) -> np.ndarray:
        return self._ad_column(0)

    @property
    def gt_alt_depths(self) -> np.ndarray:
        return self._ad_column(1)

    def _ad_column(self, column: int) -> np.ndarray:
        try:
            ad = self.format("AD")
            if ad.ndim == 1 or ad.shape[1] <= column: return np.full(len(ad), -1, dtype=np.int32)
            values = ad[:, column].astype(np.int32, copy=False)
            values[values <= -2147483640] = -1
            return values
        except KeyError: return np.full(len(self._reader.samples), -1, dtype=np.int32)

    @property
    def num_called(self) -> int: return int(np.sum(self.gt_types != UNKNOWN))
    @property
    def num_unknown(self) -> int: return int(np.sum(self.gt_types == UNKNOWN))
    @property
    def num_het(self) -> int: return int(np.sum(self.gt_types == HET))
    @property
    def num_hom_ref(self) -> int: return int(np.sum(self.gt_types == HOM_REF))
    @property
    def num_hom_alt(self) -> int: return int(np.sum(self.gt_types == HOM_ALT))
    @property
    def call_rate(self) -> float: return self.num_called / len(self._reader.samples) if self._reader.samples else 0.0
    @property
    def aaf(self) -> float:
        gt = self._reader._alleles[self._index, :, :2]
        called = gt[gt >= 0]
        return float(np.sum(called > 0) / len(called)) if len(called) else 0.0


class VCF(Iterator[Variant]):
    """Read an uncompressed VCF, with cyvcf2-compatible iteration and calls.

    ``mode`` and ``threads`` are accepted for source compatibility. Compressed VCF/BCF,
    indexed region queries, and Writer are intentionally outside this focused port.
    """

    def __init__(self, fname: str | Path, mode: str = "r", gts012: bool = False, lazy: bool = False,
                 strict_gt: bool = False, samples: list[str] | None = None, threads: int | None = None):
        if mode not in ("r", "rb"):
            raise ValueError("mojo-cyvcf2 only supports read mode")
        self.filename = str(fname)
        if self.filename.endswith((".gz", ".bgz", ".bcf")):
            raise ValueError("mojo-cyvcf2 currently supports plain-text .vcf files only")
        self.gts012, self.lazy, self.strict_gt, self.threads = gts012, lazy, strict_gt, threads
        self._data = Path(fname).read_bytes()
        self._bytes = np.frombuffer(self._data, dtype=np.uint8)
        self.raw_header, header_samples = self._read_header()
        self._file_samples = header_samples
        if samples is not None:
            missing = set(samples).difference(header_samples)
            if missing: raise ValueError(f"samples not in VCF header: {sorted(missing)}")
            self._sample_indices = [header_samples.index(s) for s in samples]
            self.samples = list(samples)
        else:
            self._sample_indices = list(range(len(header_samples)))
            self.samples = header_samples
        self._scan()
        self._cursor = 0

    def _read_header(self) -> tuple[str, list[str]]:
        text_end = self._data.find(b"#CHROM")
        if text_end < 0: raise ValueError("not a VCF: missing #CHROM header")
        line_end = self._data.find(b"\n", text_end)
        line_end = len(self._data) if line_end < 0 else line_end
        line = self._data[text_end:line_end].decode("utf-8").rstrip("\r").split("\t")
        if line[:8] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]:
            raise ValueError("invalid VCF column header")
        return self._data[:line_end + 1].decode("utf-8"), line[9:]

    def _scan(self) -> None:
        n = len(self._bytes)
        capacity = max(1, int(np.count_nonzero(self._bytes == 10)) + 1)
        starts = np.empty(capacity, dtype=np.int64)
        ends = np.empty(capacity, dtype=np.int64)
        count = lib().mcv_scan_records(addr(self._bytes), n, addr(starts), addr(ends), capacity)
        self._starts, self._ends = starts[:count], ends[:count]
        all_samples = len(self._file_samples)
        total = count * all_samples
        raw_alleles = np.empty(max(1, total * 3), dtype=np.int32)
        raw_types = np.empty(max(1, total), dtype=np.int32)
        lib().mcv_decode_gt(addr(self._bytes), addr(self._starts), addr(self._ends), count,
                            all_samples, int(self.strict_gt), addr(raw_alleles), addr(raw_types))
        decoded_alleles = raw_alleles[:total * 3].reshape(count, all_samples, 3)
        decoded_types = raw_types[:total].reshape(count, all_samples)
        if self._sample_indices == list(range(all_samples)):
            self._alleles = decoded_alleles
            self._types = decoded_types
        else:
            self._alleles = decoded_alleles[:, self._sample_indices, :]
            self._types = decoded_types[:, self._sample_indices]

    def __iter__(self) -> "VCF": self._cursor = 0; return self

    def __next__(self) -> Variant:
        if self._cursor >= len(self._starts): raise StopIteration
        record = self._variant(self._cursor)
        self._cursor += 1
        return record

    def _variant(self, index: int) -> Variant:
        fields = self._data[int(self._starts[index]):int(self._ends[index])].decode("utf-8").split("\t")
        if self._sample_indices != list(range(len(fields[9:]))):
            fields = fields[:9] + [fields[9 + i] for i in self._sample_indices]
        return Variant(self, index, tuple(fields))

    def __call__(self, region: str) -> Iterator[Variant]:
        chrom, sep, interval = region.partition(":")
        if sep:
            first, dash, last = interval.replace(",", "").partition("-")
            lo, hi = int(first), int(last) if dash else int(first)
        else:
            lo, hi = 1, 2**63 - 1
        return (v for v in (self._variant(i) for i in range(len(self._starts))) if v.CHROM == chrom and lo <= v.POS <= hi)

    def close(self) -> None: pass

    def add_info_to_header(self, header: dict[str, str]) -> None: raise NotImplementedError("header mutation is not implemented")
    def add_format_to_header(self, header: dict[str, str]) -> None: raise NotImplementedError("header mutation is not implemented")
    def add_filter_to_header(self, header: dict[str, str]) -> None: raise NotImplementedError("header mutation is not implemented")
