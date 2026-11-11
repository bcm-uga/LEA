# Text-output EOL views only; preserve sentinels, all other bytes, and goldens.
am_bytes <- function(path) {
    if (!file.exists(path)) return(NULL)
    x <- readBin(path, "raw", n = file.info(path)$size)
    if (any(x == as.raw(0))) stop("Unexpected NUL in conversion output.")
    rawToChar(x)
}
am_eol <- function(x, windows = .Platform$OS.type == "windows") {
    if (is.null(x) || !windows || !grepl("\r", x, fixed = TRUE)) return(x)
    rest <- gsub("\r\n", "", x, fixed = TRUE)
    if (grepl("\r", rest, fixed = TRUE) || grepl("\n", rest, fixed = TRUE))
        stop("Unexpected bare CR or mixed output line endings.")
    gsub("\r\n", "\n", x, fixed = TRUE)
}
am_baseline <- function(base, capture = FALSE) {
    base <- normalizePath(base)
    expected <- file.path(base, "expected")
    fixtures <- sort(list.files(file.path(base, "fixtures"), full.names = TRUE))
    root <- tempfile("ancestrymap-baseline-"); dir.create(root)
    old <- getwd(); on.exit({setwd(old); unlink(root, recursive = TRUE)})
    converters <- list(geno = LEA::ancestrymap2geno, lfmm = LEA::ancestrymap2lfmm)
    observations <- list()
    for (fixture in fixtures) {
        name <- sub("\\.ancestrymap$", "", basename(fixture))
        target <- file.path(expected, name)
        if (capture) dir.create(target, recursive = TRUE, showWarnings = FALSE)
        for (format in names(converters)) for (operation in c("default", "explicit", "native")) {
            work <- file.path(root, paste(name, format, operation, sep = "-")); dir.create(work); setwd(work)
            stopifnot(file.copy(fixture, "input.ancestrymap"))
            path <- paste0(if (operation == "explicit") "custom" else "input", ".", format)
            native <- NULL
            log <- capture.output({
                if (operation == "native") {
                    native <- .C(paste0("R_ancestrymap2", format), "input.ancestrymap", path,
                                 N = 0L, M = 0L, PACKAGE = "LEA")
                    result <- native[[2L]]
                } else if (operation == "explicit") result <- converters[[format]]("input.ancestrymap", path)
                else result <- converters[[format]]("input.ancestrymap")
            })
            golden <- file.path(target, paste0("output.", format))
            stopifnot(file.exists(path))
            if (capture && operation == "default") stopifnot(file.copy(path, golden, overwrite = TRUE))
            stopifnot(identical(am_eol(am_bytes(path)), am_eol(am_bytes(golden))))
            observations[[paste(name, format, operation, sep = "/")]] <- list(
                returned = result, stdout = log,
                counts = if (is.null(native)) NULL else c(N = native$N, M = native$M))
        }
    }
    fixture <- file.path(base, "fixtures", "matrix.ancestrymap")
    for (format in names(converters)) for (explicit in c(FALSE, TRUE))
        for (scenario in c("skip-existing", "force-existing", "no-existing")) {
            key <- paste("paths", format, explicit, scenario, sep = "/")
            work <- file.path(root, gsub("/", "-", key)); dir.create(work); setwd(work)
            stopifnot(file.copy(fixture, "input.ancestrymap"))
            path <- paste0(if (explicit) "custom" else "input", ".", format)
            if (scenario != "no-existing") writeBin(charToRaw("sentinel\n"), path)
            args <- list(input.file = "input.ancestrymap", force = scenario == "force-existing")
            if (explicit) args$output.file <- path
            log <- capture.output(result <- do.call(converters[[format]], args))
            value <- am_bytes(path)
            if (scenario != "skip-existing") value <- am_eol(value)
            observations[[key]] <- list(returned = result, stdout = log, output = value)
        }
    if (capture) dput(observations, file.path(expected, "behavior.R"))
    else stopifnot(identical(observations, dget(file.path(expected, "behavior.R"))))
    cat("ANCESTRYMAP golden compatibility:", length(fixtures) * 6L, "conversions and 12 path/force scenarios passed\n")
    invisible(observations)
}
