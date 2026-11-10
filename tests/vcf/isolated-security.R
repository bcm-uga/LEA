# Invoked ONLY by the subprocess runner, never sourced in package tests.
args <- commandArgs(TRUE)
library(LEA)
stopifnot(length(args) == 3L)
tryCatch({
    result <- if (args[1] == "geno") vcf2geno(args[2], args[3]) else vcf2lfmm(args[2])
    cat("VCF_RESULT:", result, "\n", sep = "")
}, error = function(e) {
    message <- conditionMessage(e)
    if (!nzchar(trimws(message))) quit(save = "no", status = 4L)
    cat("R_ERROR:", message, "\n", sep = "", file = stderr())
    quit(save = "no", status = 2L)
})
