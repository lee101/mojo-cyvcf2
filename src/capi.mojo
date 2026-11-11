"""Hot VCF text scanning and GT decoding exposed through a small C ABI."""

from std.sys import simd_width_of

comptime BytePtr = Pointer[UInt8, AnyOrigin[mut=True]]
comptime IntPtr = Pointer[Int, AnyOrigin[mut=True]]
comptime Int32Ptr = Pointer[Int32, AnyOrigin[mut=True]]


comptime TAB = UInt8(9)
comptime LF = UInt8(10)
comptime CR = UInt8(13)
comptime HASH = UInt8(35)
comptime COLON = UInt8(58)


def byte_at(src: BytePtr, i: Int) -> Int:
    return Int(src.unsafe_load(i))


@export("mcv_scan_records")
def mcv_scan_records(src_addr: Int, n: Int, starts_addr: Int, ends_addr: Int, capacity: Int) abi("C") -> Int:
    """Write byte ranges for non-header VCF lines and return their count."""
    var src = BytePtr(unsafe_from_address=src_addr)
    var starts = IntPtr(unsafe_from_address=starts_addr)
    var ends = IntPtr(unsafe_from_address=ends_addr)
    var line_start = 0
    var count = 0
    var i = 0
    comptime W = simd_width_of[DType.uint8]()
    while i + W <= n:
        var bytes = src.unsafe_load[width=W, alignment=1](i)
        var newlines = bytes.eq(LF)
        if newlines.reduce_or():
            comptime for lane in range(W):
                if newlines[lane]:
                    var end = i + lane
                    # A malformed line still occupies a line; the header byte is
                    # only read once we know the line is non-empty.
                    var wanted = False
                    if end > line_start:
                        wanted = src.unsafe_load(line_start) != HASH
                    if wanted:
                        if count < capacity:
                            starts.unsafe_store(count, line_start)
                            if src.unsafe_load(end - 1) == CR:
                                ends.unsafe_store(count, end - 1)
                            else:
                                ends.unsafe_store(count, end)
                        count += 1
                    line_start = end + 1
        i += W
    while i < n:
        if src.unsafe_load(i) == LF:
            if i > line_start and src.unsafe_load(line_start) != HASH:
                if count < capacity:
                    starts.unsafe_store(count, line_start)
                    if src.unsafe_load(i - 1) == CR:
                        ends.unsafe_store(count, i - 1)
                    else:
                        ends.unsafe_store(count, i)
                count += 1
            line_start = i + 1
        i += 1
    if line_start < n and src.unsafe_load(line_start) != HASH:
        if count < capacity:
            starts.unsafe_store(count, line_start)
            if src.unsafe_load(n - 1) == CR:
                ends.unsafe_store(count, n - 1)
            else:
                ends.unsafe_store(count, n)
        count += 1
    return count


def allele_type(first: Int, second: Int, missing: Bool, index: Int, strict: Int) -> Int:
    """Return the cyvcf2 genotype class for a parsed two-allele call."""
    if missing and strict != 0:
        return 2
    if missing:
        if first < 0 and second < 0:
            return 2
        if first == 0 or second == 0:
            return 0
        return 1
    if index == 0:
        return 0 if first == 0 else 3
    if first < 0 or second < 0:
        return 2
    if first == 0 and second == 0:
        return 0
    if first != second:
        return 1
    return 3


def parse_allele(src: BytePtr, begin: Int, finish: Int, strict: Int, alleles: Int32Ptr, dst: Int) -> Int:
    comptime SLASH = UInt8(47)
    comptime PIPE = UInt8(124)
    comptime DOT = UInt8(46)
    var first = -1
    var second = -1
    var phased = 0
    var length = finish - begin
    if length == 1:
        # Haploid call: one allele, or a lone missing value.
        var c0 = src.unsafe_load(begin)
        if c0 == DOT:
            alleles.unsafe_store(dst, Int32(-1))
            alleles.unsafe_store(dst + 1, Int32(-1))
            alleles.unsafe_store(dst + 2, Int32(0))
            return allele_type(-1, -1, True, 0, strict)
        if c0 >= 48 and c0 <= 57:
            first = Int(c0) - 48
            alleles.unsafe_store(dst, Int32(first))
            alleles.unsafe_store(dst + 1, Int32(-1))
            alleles.unsafe_store(dst + 2, Int32(0))
            return 0 if first == 0 else 3
    elif length == 3:
        # Diploid, single-digit call: ``a/b``, ``a|b`` or a missing allele.  The
        # fixed length proves both allele bytes are in range.
        var separator = src.unsafe_load(begin + 1)
        if separator == SLASH or separator == PIPE:
            var c0 = src.unsafe_load(begin)
            var c2 = src.unsafe_load(begin + 2)
            var left = (c0 >= 48 and c0 <= 57) or c0 == DOT
            var right = (c2 >= 48 and c2 <= 57) or c2 == DOT
            if left and right:
                phased = 1 if separator == PIPE else 0
                if c0 != DOT:
                    first = Int(c0) - 48
                if c2 != DOT:
                    second = Int(c2) - 48
                alleles.unsafe_store(dst, Int32(first))
                alleles.unsafe_store(dst + 1, Int32(second))
                alleles.unsafe_store(dst + 2, Int32(phased))
                var missing = first < 0 or second < 0
                return allele_type(first, second, missing, 1, strict)

    # General path: multi-digit alleles, polyploid calls, and malformed text.
    var value = 0
    var seen_digit = False
    var index = 0
    first = -1
    second = -1
    phased = 0
    var missing_field = False
    var i = begin
    while i < finish:
        var c = byte_at(src, i)
        if c >= 48 and c <= 57:
            value = value * 10 + c - 48
            seen_digit = True
        elif c == Int(SLASH) or c == Int(PIPE):
            if c == Int(PIPE):
                phased = 1
            if not seen_digit:
                missing_field = True
            elif index == 0:
                first = value
            else:
                second = value
            index += 1
            value = 0
            seen_digit = False
        else:
            missing_field = True
        i += 1
    if not seen_digit:
        missing_field = True
    elif index == 0:
        first = value
    else:
        second = value
    if index == 0:
        second = -1
    alleles.unsafe_store(dst, Int32(first))
    alleles.unsafe_store(dst + 1, Int32(second))
    alleles.unsafe_store(dst + 2, Int32(phased))
    return allele_type(first, second, missing_field, index, strict)



def decode_record(src: BytePtr, starts: IntPtr, ends: IntPtr, r: Int, samples: Int, strict: Int, alleles: Int32Ptr, types: Int32Ptr):
    var finish = ends.unsafe_load(r)
    var field_begin = starts.unsafe_load(r)
    var format_begin = finish
    var format_end = finish
    var sample_begin = finish
    # Walk the eight fixed leading columns with a tab counter, then treat the
    # ninth column as FORMAT.  A record without a ninth tab has no samples and
    # its FORMAT runs to `finish`.
    var column = 0
    var lead = field_begin
    while lead < finish:
        if src.unsafe_load(lead) == TAB:
            column += 1
            field_begin = lead + 1
            if column == 8:
                break
        lead += 1
    if column == 8 and field_begin < finish:
        format_begin = field_begin
        var ninth = field_begin
        while ninth < finish and src.unsafe_load(ninth) != TAB:
            ninth += 1
        format_end = ninth
        if ninth < finish:
            sample_begin = ninth + 1

    var gt_column = -1
    if format_begin < format_end:
        var key_begin = format_begin
        var key = 0
        var key_index = format_begin
        while key_index <= format_end:
            if key_index == format_end or src.unsafe_load(key_index) == COLON:
                if key_index - key_begin == 2 and src.unsafe_load(key_begin) == 71 and src.unsafe_load(key_begin + 1) == 84:
                    gt_column = key
                    break
                key += 1
                key_begin = key_index + 1
            key_index += 1

    var sample = 0
    var type_slot = r * samples
    while sample < samples:
        var sample_end = sample_begin
        while sample_end < finish and src.unsafe_load(sample_end) != TAB:
            sample_end += 1
        var gt_begin = sample_end
        var gt_end = sample_end
        if gt_column == 0:
            # GT is the first FORMAT value, so the call starts the sample field.
            gt_begin = sample_begin
            gt_end = gt_begin
            while gt_end < sample_end and src.unsafe_load(gt_end) != COLON:
                gt_end += 1
        elif gt_column > 0:
            var value_begin = sample_begin
            var value_column = 0
            var cursor = sample_begin
            while cursor < sample_end and value_column < gt_column:
                if src.unsafe_load(cursor) == COLON:
                    value_column += 1
                    value_begin = cursor + 1
                cursor += 1
            if value_column == gt_column:
                gt_begin = value_begin
                gt_end = gt_begin
                while gt_end < sample_end and src.unsafe_load(gt_end) != COLON:
                    gt_end += 1
        var allele_slot = type_slot * 3
        types.unsafe_store(type_slot, Int32(parse_allele(src, gt_begin, gt_end, strict, alleles, allele_slot)))
        type_slot += 1
        sample_begin = sample_end + 1
        sample += 1


@export("mcv_decode_gt")
def mcv_decode_gt(src_addr: Int, starts_addr: Int, ends_addr: Int, records: Int, samples: Int, strict: Int, alleles_addr: Int, types_addr: Int) abi("C"):
    """Decode GT values into (allele1, allele2, phased) and cyvcf2 types."""
    var src = BytePtr(unsafe_from_address=src_addr)
    var starts = IntPtr(unsafe_from_address=starts_addr)
    var ends = IntPtr(unsafe_from_address=ends_addr)
    var alleles = Int32Ptr(unsafe_from_address=alleles_addr)
    var types = Int32Ptr(unsafe_from_address=types_addr)
    var r = 0
    while r < records:
        decode_record(src, starts, ends, r, samples, strict, alleles, types)
        r += 1
