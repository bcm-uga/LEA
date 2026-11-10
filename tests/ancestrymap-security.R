library(LEA)

local({
    root <- tempfile("ancestrymap-security-")
    dir.create(root)
    on.exit(unlink(root, recursive = TRUE))
    converters <- list(geno = ancestrymap2geno, lfmm = ancestrymap2lfmm)
    checks <- 0L
    run <- function(records, expected = NULL, pattern = NULL, newline = TRUE) {
        input <- file.path(root, "input.ancestrymap")
        bytes <- paste(records, collapse = "\n")
        if (newline && length(records)) bytes <- paste0(bytes, "\n")
        writeBin(charToRaw(bytes), input)
        for (format in names(converters)) {
            output <- file.path(root, paste0("output.", format))
            # Invalid input must not overwrite an existing output.
            writeLines("sentinel", output)
            result <- tryCatch(converters[[format]](input, output, force = TRUE),
                               error = identity)
            if (is.null(pattern)) {
                stopifnot(!inherits(result, "error"), identical(result, output))
                observed <- readChar(output, file.info(output)$size, useBytes = TRUE)
                stopifnot(identical(observed, expected[[format]]))
            } else {
                stopifnot(inherits(result, "error"),
                          grepl(pattern, conditionMessage(result)),
                          identical(readLines(output), "sentinel"))
            }
            checks <<- checks + 1L
        }
    }
    multiple <- list(geno = "01\n29\n", lfmm = "0 2\n1 9\n")
    single <- list(geno = "01\n", lfmm = "0\n1\n")
    oversized <- "SNP name too long.*maximum is 511 bytes"
    for (n in c(1L, 510L, 511L, 512L, 513L, 600L)) {
        name <- strrep("a", n)
        for (first in c(TRUE, FALSE)) {
            names <- if (first) c(name, "b") else c("b", name)
            records <- c(paste(names[1], "ind1 0"), paste(names[1], "ind2 1"),
                         paste(names[2], "ind1 2"), paste(names[2], "ind2 9"))
            run(records, multiple, if (n >= 512L) oversized else NULL)
        }
        # The third locus is beyond the counting pass and exercises the reader.
        run(c("b ind1 0", "b ind2 1", "c ind1 2", "c ind2 9",
              paste(name, "ind1 0"), paste(name, "ind2 1")),
            list(geno = "01\n29\n01\n", lfmm = "0 2 0\n1 9 1\n"),
            if (n >= 512L) oversized else NULL)
        if (n <= 511L) {
            for (newline in c(TRUE, FALSE))
                run(c(paste(name, "ind1 0"), paste(name, "ind2 1")),
                    single, newline = newline)
        }
    }
    # Distinct names differing only in byte 511 must not merge.
    prefix <- strrep("a", 510L)
    run(c(paste0(prefix, "x ind1 0"), paste0(prefix, "x ind2 1"),
          paste0(prefix, "y ind1 2"), paste0(prefix, "y ind2 9")), multiple)
    run(c("a ind1 0", "a ind2 1", "b ind1 2", "b ind2 9"),
        multiple, newline = FALSE)
    run("a ind1 2", list(geno = "2\n", lfmm = "2\n"), newline = FALSE)
    run(c("a\tind1\t0\r", "a\tind2\t1\r"), single)
    run(character(), pattern = "Empty ANCESTRYMAP")
    for (records in list("", " \t", "a", "a ind1", c("a ind1 0", ""),
                         c("a ind1 0", "b ind1")))
        run(records, pattern = "Invalid ANCESTRYMAP record")
    # An incomplete final locus and unequal group sizes must be controlled errors.
    run(c("a i 0", "a j 1", "b i 2"), pattern = "incorrect number of lines")
    run(c("a i 0", "a j 1", "b i 2", "b j 9", "b k 0", "c i 1"),
        pattern = "Invalid ANCESTRYMAP group sizes")
    # Bundled fixture has a descriptive header; convert its data records and
    # compare against an independent matrix conversion.
    example <- read.table(gzfile(system.file("data", "example_ancestrymap.txt.gz",
                                            package = "LEA")),
                          colClasses = "character", skip = 1L)
    n <- sum(example[[1]] == example[[1]][1])
    genotypes <- matrix(as.integer(example[[3]]), nrow = n)
    expected <- list(
        geno = paste0(paste(apply(genotypes, 2, paste0, collapse = ""),
                            collapse = "\n"), "\n"),
        lfmm = paste0(paste(apply(genotypes, 1, paste, collapse = " "),
                            collapse = "\n"), "\n"))
    run(apply(example, 1, paste, collapse = " "), expected)
    cat("ANCESTRYMAP security regression checks passed:", checks, "\n")
})
