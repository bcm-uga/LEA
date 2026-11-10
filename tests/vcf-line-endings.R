# Exercise the Windows comparison branch on every host without claiming a
# Windows runtime test. Stored golden files are never rewritten.
base <- if (dir.exists("vcf/fixtures")) "vcf" else "tests/vcf"
source(file.path(base, "baseline-helpers.R"))
canonical <- "1 A \n2\tB \n"
stopifnot(identical(vcf_normalize_output_eol(canonical, TRUE), canonical),
          identical(vcf_normalize_output_eol("1 A \r\n2\tB \r\n", TRUE), canonical),
          identical(vcf_normalize_output_eol("", TRUE), ""),
          is.null(vcf_normalize_output_eol(NULL, TRUE)))
for (invalid in c("a\rb\n", "a\r\r\n", "a\r\nb\n", "a\nb\r\n")) {
    failure <- tryCatch(vcf_normalize_output_eol(invalid, TRUE), error = identity)
    stopifnot(inherits(failure, "error"))
}
# Genuine format/content changes must still compare unequal.
for (changed in c("1 A\r\n2\tB \r\n", "1 A \r\n2 B \r\n", "2\tB \r\n1 A \r\n",
                  "1 A \r\n2\tB ", "1 A \r\n2\tB \r\n\r\n", "0 A \r\n2\tB \r\n"))
    stopifnot(!identical(vcf_normalize_output_eol(changed, TRUE), canonical))
stopifnot(!identical(vcf_normalize_output_eol("\r\n", TRUE), ""),
          identical(vcf_normalize_output_eol("a\r\n", FALSE), "a\r\n"),
          !identical(vcf_normalize_output_eol("a\r\n", FALSE), "a\n"))
# The sentinel is intentionally binary on every platform, including Windows.
path <- tempfile("vcf-sentinel-")
writeBin(charToRaw("sentinel\n"), path)
stopifnot(identical(readBin(path, "raw", n = 100), charToRaw("sentinel\n")))
stopifnot(identical(vcf_read_output_bytes(path), "sentinel\n"))
writeBin(c(charToRaw("sentinel\n"), as.raw(0), charToRaw("unexpected")), path)
failure <- tryCatch(vcf_read_output_bytes(path), error = identity)
stopifnot(inherits(failure, "error"))
# A CRLF checkout-style expectation is normalized only in memory.
writeBin(charToRaw("1 A \r\n2\tB \r\n"), path)
before <- readBin(path, "raw", n = 100)
stopifnot(identical(vcf_normalize_output_eol(vcf_read_output_bytes(path), TRUE), canonical),
          identical(readBin(path, "raw", n = 100), before))
unlink(path)
cat("VCF line-ending checks passed (Windows branch simulated; not Windows execution)\n")
