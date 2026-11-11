library(LEA)
local({
    root <- tempfile("vcf-controls-"); dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    input <- file.path(root, "input.vcf")
    output <- file.path(root, "output.geno")
    header <- charToRaw("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
    tail <- rep(charToRaw("A"), 1024 * 1024)
    forbidden <- c(0:8, 11:12, 14:31, 127)
    for (byte in forbidden) {
        writeBin(c(header, as.raw(byte), tail), input)
        failure <- tryCatch(vcf2geno(input, output), error = identity)
        stopifnot(inherits(failure, "error"), grepl("VCF error at line 2", conditionMessage(failure)))
    }
    cat(length(forbidden), "real-R control-byte/long-tail rejection checks passed\n")
})
