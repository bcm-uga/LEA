library(LEA)
base <- if (dir.exists("vcf/fixtures")) "vcf" else "tests/vcf"
source(file.path(base, "baseline-helpers.R"))
vcf_baseline(base)
