# ANCESTRYMAP compatibility decisions (recorded before implementation)

Captured from unmodified c56f629, LEA 3.25.1. The isolated baseline library was
built from that commit during the cleanup audit; source equality was verified
before capture. Do not regenerate the expected files after implementing fixes.

Run `R_LIBS=/path/to/library Rscript tests/ancestrymap-compatibility.R`.
Nine fixtures cover 54 conversions (default/explicit paths and native counts for
both formats), plus 12 existing/missing output and force scenarios. Returned paths,
stdout including warnings/counts, and exact matrix outputs are frozen. Binary
sentinels and scoped Windows EOL comparison views preserve all other bytes.
The bundled fixture removes only its descriptive header, as existing tests do.

Policy chosen for the cleanup:

* Keep ignoring fields after the third as historical annotations. This is safe
  when the entire physical record has been read/checked; do not silently drop
  separate physical records or accept embedded NUL/control bytes.
* Keep representable decimal integers outside 0/1/2/9, with the existing once-per-
  conversion warning and exact output spelling after conversion to int. Preserve
  leading signs/zero padding. Reject junk, fractional tokens, and int/long overflow
  rather than retaining atoi's silent mis-conversions/undefined behavior.
* Preserve existing warning text, ordering in captured stdout, and historical
  warning coordinates (even their legacy i/j convention). Defer printing until
  resources are released; do not change valid conversion outputs to fix cosmetic
  diagnostics. Parsing errors may use reliable physical-line context.
* Keep the 511-byte SNP-name bound, sample/variant order, ignored ID semantics,
  GENO/LFMM orientation, separators and platform text-mode output behavior.
* Keep output opening after input validation; invalid data must not overwrite
  existing output. Preserve documented default/explicit path and force behavior.
* Reject invalid filename vectors; length-zero/NA/multiple filenames are not
  valid uses of the documented single-string interface.
* Explicitly reject counts beyond public int/INT_MAX limits rather than extending
  multi-billion-record support. Use checked size_t arithmetic for storage/indexing.

The manifest pins every fixture/golden SHA-256. Existing 66 native and 76 R tests
are separate, unchanged regression suites. No post-fix golden capture is permitted.
