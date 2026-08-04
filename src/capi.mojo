"""Hot VCF text scanning and GT decoding exposed through a small C ABI."""

from std.sys import simd_width_of

comptime BytePtr = UnsafePointer[UInt8, AnyOrigin[mut=True]]
comptime IntPtr = UnsafePointer[Int, AnyOrigin[mut=True]]
comptime Int32Ptr = UnsafePointer[Int32, AnyOrigin[mut=True]]


def byte_at(src: BytePtr, i: Int) -> Int:
    return Int(src.load(i))


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
        var bytes = src.load[width=W](i)
        comptime for lane in range(W):
            if bytes[lane] == 10:
                var end = i + lane
                if end > line_start and byte_at(src, line_start) != 35:
                    if count < capacity:
                        starts.store(count, line_start)
                        if end > line_start and byte_at(src, end - 1) == 13:
                            ends.store(count, end - 1)
                        else:
                            ends.store(count, end)
                    count += 1
                line_start = end + 1
        i += W
    while i < n:
        if byte_at(src, i) == 10:
            if i > line_start and byte_at(src, line_start) != 35:
                if count < capacity:
                    starts.store(count, line_start)
                    if i > line_start and byte_at(src, i - 1) == 13:
                        ends.store(count, i - 1)
                    else:
                        ends.store(count, i)
                count += 1
            line_start = i + 1
        i += 1
    if line_start < n and byte_at(src, line_start) != 35:
        if count < capacity:
            starts.store(count, line_start)
            if n > line_start and byte_at(src, n - 1) == 13:
                ends.store(count, n - 1)
            else:
                ends.store(count, n)
        count += 1
    return count


def parse_allele(src: BytePtr, begin: Int, finish: Int, strict: Int, alleles: Int32Ptr, dst: Int) -> Int:
    var first = -1
    var second = -1
    var value = 0
    var seen_digit = False
    var index = 0
    var phase = 0
    var missing = False
    var i = begin
    while i < finish:
        var c = byte_at(src, i)
        if c >= 48 and c <= 57:
            value = value * 10 + c - 48
            seen_digit = True
        elif c == 47 or c == 124:
            if c == 124:
                phase = 1
            if not seen_digit:
                missing = True
            elif index == 0:
                first = value
            else:
                second = value
            index += 1
            value = 0
            seen_digit = False
        else:
            missing = True
        i += 1
    if not seen_digit:
        missing = True
    elif index == 0:
        first = value
    else:
        second = value
    if index == 0:
        second = -1
    alleles.store(dst, Int32(first))
    alleles.store(dst + 1, Int32(second))
    alleles.store(dst + 2, Int32(phase))
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


def decode_record(src: BytePtr, starts: IntPtr, ends: IntPtr, r: Int, samples: Int, strict: Int, alleles: Int32Ptr, types: Int32Ptr):
    var finish = ends.load(r)
    var field_begin = starts.load(r)
    var format_begin = finish
    var format_end = finish
    var sample_begin = finish
    var column = 0
    var i = field_begin
    while i < finish:
        if byte_at(src, i) == 9:
            if column == 7:
                # FORMAT is allowed to be the final column when a VCF has no
                # samples.  In that case format_end remains `finish`.
                format_begin = i + 1
            elif column == 8:
                format_begin = field_begin
                format_end = i
                sample_begin = i + 1
                break
            column += 1
            field_begin = i + 1
        i += 1

    var gt_column = -1
    if format_begin < format_end:
        var key_begin = format_begin
        var key = 0
        i = format_begin
        while i <= format_end:
            if i == format_end or byte_at(src, i) == 58:
                if i - key_begin == 2 and byte_at(src, key_begin) == 71 and byte_at(src, key_begin + 1) == 84:
                    gt_column = key
                    break
                key += 1
                key_begin = i + 1
            i += 1

    var sample = 0
    while sample < samples:
        var sample_end = sample_begin
        while sample_end < finish and byte_at(src, sample_end) != 9:
            sample_end += 1
        var gt_begin = sample_end
        var gt_end = sample_end
        if gt_column >= 0:
            var value_begin = sample_begin
            var value_column = 0
            i = sample_begin
            while i < sample_end and value_column < gt_column:
                if byte_at(src, i) == 58:
                    value_column += 1
                    value_begin = i + 1
                i += 1
            if value_column == gt_column:
                gt_begin = value_begin
                gt_end = gt_begin
                while gt_end < sample_end and byte_at(src, gt_end) != 58:
                    gt_end += 1
        var slot = (r * samples + sample) * 3
        types.store(r * samples + sample, Int32(parse_allele(src, gt_begin, gt_end, strict, alleles, slot)))
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
