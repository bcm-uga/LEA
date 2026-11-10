# ANCESTRYMAP security fix report

Date: 2026-11-10. Package version: **3.25.0**, unchanged.
Branch: `fix/ancestrymap-buffer-overflow`.
Base commit: `89a169585b63bf4b1a6672262f5ebe31e9a0601d`.

## Changes

Only `src/convert/ancestrymap.c`, `src/convert/ancestrymap.h`, new ANCESTRYMAP regression tests, and this report were changed/added for implementation. The earlier audit remains unchanged. No R source, DESCRIPTION, NAMESPACE, VCF conversion, LFMM2 implementation, benchmarks, or other production code was modified. Nothing was committed or pushed.

* Defined `ANCESTRYMAP_SNP_NAME_CAPACITY` as 512 bytes, including the terminating NUL, and used it for all ANCESTRYMAP name arrays. The parser's destination-buffer contract is documented in its header.
* Added a shared internal validator, called before copying in both `read_line_ancestrymap` and `nb_ind_ancestrymap`. NULL tokens raise a descriptive error. Tokens of 512 bytes or more raise an R error identifying the file, parsing position, and 511-byte limit. Valid names are copied exactly; no truncation is performed.
* Reworked the first-locus counting loop to tokenize only successful `fgets` results. It counts each first-locus record once, stops at the first different locus or EOF, checks stream errors, and rejects empty files. Valid single-locus files work with or without a final newline.
* Missing individual/genotype columns now produce descriptive R errors. Genotype conversion and warning behavior remain unchanged for existing inputs.
* Added checks before matrix writes to reject inconsistent individual/locus groups. This is limited to ANCESTRYMAP and prevents malformed group sizes from indexing outside the allocated data matrix.
* Replaced the ANCESTRYMAP dimension/group error paths with descriptive non-NULL R error messages. Regression testing demonstrated that the old `error(NULL)` path for an incomplete locus crashes R 4.6.1; keeping that path would not satisfy controlled-error handling. Shared I/O helpers were not changed.

## Tests added

`tests/ancestrymap-security.R` uses base R and the installed LEA package, with no new package dependencies. It makes **76 conversion assertions** across both `ancestrymap2geno` and `ancestrymap2lfmm`:

* SNP-name lengths 1, 510, 511, 512, 513, and 600 bytes; first, second, and third locus placement. The third locus reaches the main reader after counting has stopped.
* Single-locus and multiple-locus inputs, one individual, optional final newline, tabs, and CRLF.
* Distinct 511-byte names that differ only in their final byte, to detect accidental truncation/merging.
* Empty/blank inputs, missing columns, EOF in an incomplete locus, and inconsistent group sizes.
* Exact GENO/LFMM bytes for known matrices and the bundled example's data records. The bundled text fixture has a descriptive header; the test excludes that header because the native format reader expects data records.
* Rejected inputs leave an existing output file's sentinel content unchanged.

`tests/native/ancestrymap-security.c` and `tests/native/run-ancestrymap-asan.py` provide a repeatable standalone sanitizer test. The driver links the production ANCESTRYMAP parser, GENO writer, numeric I/O, and shared file/error helpers unchanged. Only R logging and R error exit are substituted. It executes each case in a separate process and requires both the intended exit status/message and absence of sanitizer diagnostics. It checks direct line parsing as well as both complete native conversions, including malformed group sizes and EOF. All build/input/output artifacts are temporary.

## Verification results

Environment: macOS arm64; Apple clang 21.0.0; R 4.6.1.

| Verification | Result |
| --- | --- |
| New native ASan suite | **66 checks passed**, no ASan diagnostics. |
| Original audit reproducer, rebuilt against corrected source | Both `line` and `count` modes: 511 bytes succeed; 512 and 600 bytes exit through the controlled error path (status 2), with no ASan diagnostics. Six cases passed. |
| Fresh package installation | Passed in an isolated temporary package copy/library; all native code rebuilt initially with `R CMD INSTALL --preclean`. |
| R-level conversion suite | **76 checks passed** against the newly installed package. |
| `git diff --check` | Passed. |
| `R CMD check --no-manual --no-build-vignettes` | Completed, exit status 0; 0 errors, 2 documentation warnings. All tests, examples, and vignette R code passed. |

The original audit reproducer's error stub exits without formatting R error text; the new native suite and actual R-level tests additionally verify the informative messages.

Commands, from the working repository unless noted:

```sh
python3 tests/native/run-ancestrymap-asan.py
R_LIBS=/private/tmp/lea-security-fix-5qydho9j/library Rscript tests/ancestrymap-security.R
```

The ASan suite compiles with `-g -O1 -std=c11 -fsanitize=address -fno-omit-frame-pointer`. Leak detection is disabled for the standalone driver because its R-error substitute exits the process rather than unwinding through R. All address-access checks remain active.

Build/install/check work was performed under `/private/tmp/lea-security-fix-5qydho9j/`, leaving the working tree's compiled artifacts untouched. Logs there include `install.log`, `regression.log`, `build.log`, `check.log`, and `LEA.Rcheck/00check.log`. Temporary artifacts are not committed and may later disappear. The original reproducer's fixed-run logs are under `/private/tmp/lea-ancestrymap-audit/fixed-*.log`.

## Package check details

In the temporary build directory, ran:

```sh
R CMD build --no-build-vignettes LEA
_R_CHECK_CRAN_INCOMING_REMOTE_=false R CMD check --no-manual --no-build-vignettes LEA_3.25.0.tar.gz
```

The check ran all 76 new R conversion checks successfully. Its two warnings concern documentation artifacts: files in `vignettes/` without built files in `inst/doc/`, and `LEA.Rnw` without a corresponding built PDF/HTML. Vignette rebuilding was deliberately skipped; these warnings do not identify a parser/test failure. Dependency-index lookups also printed unavailable-repository messages, but the dependency check passed using installed dependencies. No production or metadata changes were made to address documentation/environment issues.

## Remaining limitations

* The independent VCF metadata/path issues identified in the audit remain untouched, as requested. This patch is not a general conversion-parser hardening effort.
* Existing fixed line-buffer sizes, dimension/allocation arithmetic, genotype numeric parsing, and extra-column behavior were not redesigned. Arbitrarily long physical records and extremely large data sets retain their existing limitations.
* As in the original implementation, errors raised during parsing can bypass `fclose`/`free` for an open input or allocated matrix. General R unwind-safe resource management is outside this minimal patch. ASan results establish absence of invalid memory accesses in the exercised cases, not absence of leaks.
* File stream error branches were checked in source but not exercised with injected I/O failures. ASan tests use actual native conversion code with an R error substitute; the R-level suite tests real R error propagation, but a sanitizer-enabled R process was not built or tested.
* Check settings omit PDF manual generation and vignette rebuilding; those documentation toolchains are outside this security correction.

Implementation is ready for review; no commit or push has been made.
