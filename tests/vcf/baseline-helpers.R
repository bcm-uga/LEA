# Canonicalize text-output comparison views only on Windows. Reject mixed or
# bare-CR endings; preserve every other byte and never rewrite golden files.
vcf_normalize_output_eol <- function(bytes, windows = .Platform$OS.type == "windows") {
    if (is.null(bytes) || !windows || !grepl("\r", bytes, fixed = TRUE)) return(bytes)
    remainder <- gsub("\r\n", "", bytes, fixed = TRUE)
    if (grepl("\r", remainder, fixed = TRUE))
        stop("Unexpected bare CR in generated VCF output.")
    if (grepl("\n", remainder, fixed = TRUE))
        stop("Mixed LF and CRLF in generated VCF output.")
    gsub("\r\n", "\n", bytes, fixed = TRUE)
}

# Binary reads must not silently truncate unexpected NUL bytes as readChar can.
vcf_read_output_bytes <- function(path) {
    if (!file.exists(path)) return(NULL)
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    if (any(bytes == as.raw(0))) stop("Unexpected NUL in VCF compatibility output.")
    rawToChar(bytes)
}

# Shared by the read-only regression test and the explicit capture command.
vcf_baseline <- function(base, capture = FALSE) {
    base <- normalizePath(base)
    fixtures <- sort(list.files(file.path(base, "fixtures"), "\\.vcf$", full.names = TRUE))
    expected <- file.path(base, "expected")
    root <- tempfile("vcf-compatibility-")
    dir.create(root)
    old <- getwd()
    on.exit({setwd(old); unlink(root, recursive = TRUE)}, add = TRUE)
    read_bytes <- vcf_read_output_bytes
    read_generated <- function(path) vcf_normalize_output_eol(read_bytes(path))
    operations <- c("geno-default", "geno-explicit", "lfmm-default", "lfmm-explicit", "native")
    observations <- list()
    for (fixture in fixtures) {
        name <- sub("\\.vcf$", "", basename(fixture))
        target <- file.path(expected, name)
        if (capture) dir.create(target, recursive = TRUE, showWarnings = FALSE)
        for (operation in operations) {
            work <- file.path(root, paste(name, operation, sep = "-"))
            dir.create(work); setwd(work)
            stopifnot(file.copy(fixture, "input.vcf"))
            native <- NULL
            log <- capture.output({
                if (operation == "native") {
                    native <- .C("R_vcf2geno", as.character("input.vcf"),
                                 as.character("input.geno"), N = as.integer(0),
                                 M = as.integer(0), PACKAGE = "LEA")
                    returned <- native[[2L]]
                } else if (operation == "geno-default") returned <- LEA::vcf2geno("input.vcf")
                else if (operation == "geno-explicit") returned <- LEA::vcf2geno("input.vcf", "custom.geno")
                else if (operation == "lfmm-default") returned <- LEA::vcf2lfmm("input.vcf")
                else returned <- LEA::vcf2lfmm("input.vcf", "custom.lfmm")
            })
            stem <- if (operation == "geno-explicit") "custom" else "input"
            suffixes <- c("geno", "vcfsnp", "removed")
            if (grepl("^lfmm", operation)) suffixes <- c(suffixes, "lfmm")
            for (suffix in suffixes) {
                actual <- paste0(stem, ".", suffix)
                golden <- file.path(target, paste0("output.", suffix))
                stopifnot(file.exists(actual))
                if (capture && (operation == "geno-default" ||
                                (operation == "lfmm-default" && suffix == "lfmm")))
                    stopifnot(file.copy(actual, golden, overwrite = TRUE))
                # A Windows checkout can translate golden-file line endings too;
                # compare canonical views, leaving stored golden bytes untouched.
                stopifnot(identical(read_generated(actual), read_generated(golden)))
            }
            if (operation == "lfmm-explicit") stopifnot(!file.exists("custom.lfmm"))
            observations[[paste(name, operation, sep = "/")]] <- list(
                returned = returned, stdout = log,
                counts = if (is.null(native)) NULL else c(N = native$N, M = native$M))
        }
    }
    # Record wrapper behavior independently of input-content assertions.
    fixture <- fixtures[basename(fixtures) == "diploid-tabs.vcf"]
    scenarios <- list(
        geno_default_skip = list(fun = "geno", output = NULL, force = FALSE, existing = "input.geno", preserve = TRUE),
        geno_explicit_skip = list(fun = "geno", output = "custom.geno", force = FALSE, existing = "custom.geno", preserve = TRUE),
        geno_force = list(fun = "geno", output = NULL, force = TRUE, existing = "input.geno"),
        lfmm_default_geno_skip = list(fun = "lfmm", output = NULL, force = FALSE, existing = "input.geno", preserve = TRUE),
        lfmm_default_lfmm_exists = list(fun = "lfmm", output = NULL, force = FALSE, existing = "input.lfmm"),
        lfmm_explicit_skip = list(fun = "lfmm", output = "custom.lfmm", force = FALSE, existing = "custom.lfmm", preserve = TRUE),
        lfmm_explicit_force = list(fun = "lfmm", output = "custom.lfmm", force = TRUE, existing = "custom.lfmm", preserve = TRUE),
        lfmm_explicit_missing = list(fun = "lfmm", output = "custom.lfmm", force = FALSE, existing = "input.geno"))
    for (name in names(scenarios)) {
        s <- scenarios[[name]]
        work <- file.path(root, name); dir.create(work); setwd(work)
        stopifnot(file.copy(fixture, "input.vcf"))
        writeBin(charToRaw("sentinel\n"), s$existing)
        args <- list(input.file = "input.vcf", force = s$force)
        if (!is.null(s$output)) args$output.file <- s$output
        log <- capture.output(returned <- do.call(if (s$fun == "geno") LEA::vcf2geno else LEA::vcf2lfmm, args))
        files <- sort(setdiff(list.files(), "input.vcf"))
        # Skipped/preserved sentinels are binary test data, not generated text.
        read_snapshot <- function(file) {
            if (isTRUE(s$preserve) && identical(file, s$existing)) read_bytes(file)
            else read_generated(file)
        }
        observations[[paste0("paths/", name)]] <- list(returned = returned, stdout = log,
                                                      files = setNames(lapply(files, read_snapshot), files))
    }
    if (capture) dput(observations, file.path(expected, "behavior.R"))
    else stopifnot(identical(observations, dget(file.path(expected, "behavior.R"))))
    cat("VCF compatibility baseline passed:", length(fixtures), "fixtures,",
        length(fixtures) * length(operations), "conversion runs,",
        length(scenarios), "path/force scenarios\n")
    invisible(observations)
}
