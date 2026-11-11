# Explicit maintenance command; never sourced by normal package tests.
# Run only against the unmodified c56f629 baseline; do not regenerate after fixes.
# R_LIBS=/path/to/baseline/library Rscript tests/vcf/capture-baseline.R
library(LEA)
stopifnot(as.character(packageVersion("LEA")) == "3.25.1")
expected.source <- dget("tests/vcf/baseline-source-md5.R")
stopifnot(identical(tools::md5sum(names(expected.source)), expected.source))
source("tests/vcf/baseline-helpers.R")
vcf_baseline("tests/vcf", capture = TRUE)
