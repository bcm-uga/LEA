library(LEA)
base <- if (dir.exists("ancestrymap/fixtures")) "ancestrymap" else "tests/ancestrymap"
source(file.path(base, "baseline-helpers.R"))
am_baseline(base)
# Source-level Windows comparison tests, not a Windows runtime claim.
stopifnot(identical(am_eol("0 1 \r\n9\t2\r\n", TRUE), "0 1 \n9\t2\n"),
          !identical(am_eol("0 1\r\n9\t2\r\n", TRUE), "0 1 \n9\t2\n"),
          identical(am_eol("", TRUE), ""),
          identical(am_eol("0\r\n", FALSE), "0\r\n"))
for (bad in c("0\r1\n","0\r\n1\n"))
    stopifnot(inherits(tryCatch(am_eol(bad,TRUE),error=identity),"error"))
