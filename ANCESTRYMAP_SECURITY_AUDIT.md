# ANCESTRYMAP security audit — LEA 3.25.0

Date: 2026-11-10. Branch: `fix/ancestrymap-buffer-overflow`.
Audited commit: `89a169585b63bf4b1a6672262f5ebe31e9a0601d`.
Reporter: Hong Qin, Old Dominion University (attribution supplied with the report).

## Conclusion and scope

The reported stack-buffer overflow is present in the checked-out 3.25.0 source. A 600-byte SNP name causes `strcpy` to copy 601 bytes, including the terminating NUL, into a 512-byte destination. There is also an earlier unchecked copy in `nb_ind_ancestrymap`; fixing only the reported line would leave the public conversion functions vulnerable. Both ANCESTRYMAP output formats are affected. VCF conversion has analogous heap-buffer overflows and a separate output-path stack-buffer risk.

This audit inspected all C/header files in `src/convert/`, the relevant R/native wrappers, registration, and shared I/O helpers. Temporary AddressSanitizer harnesses compiled the original parser sources without editing them. No C/R source, package metadata, existing tests, LFMM2 implementation, or benchmarks were modified. No commit or push was performed. Historical 3.23.0 sources and the researcher's original reproducer were not independently examined; the current source and independent reproductions establish the defect in 3.25.0.

## 1. Exact cause and name limit

In `src/convert/ancestrymap.c:97–145`, `read_ancestrymap` allocates stack arrays `tmp[512]`, `ref[512]`, and `szbuff[1000]`. It passes `tmp` as the `name` argument to `read_line_ancestrymap` at line 119. That function takes the first whitespace-delimited token using `strtok(szbuff, SEP)` and executes `strcpy(name, token)` at line 160, without checking its length or receiving the destination capacity.

`SEP` is `" \t\r\n"` (`src/io/read.h:14`). The `fgets` bound protects the 1000-byte line buffer, not the smaller name buffer. A 600-byte first token plus ordinary individual/genotype fields fits within that line buffer. The copy writes 89 bytes beyond a 512-byte destination (601 minus 512). A name of exactly 512 bytes already causes a one-byte NUL overflow.

**The maximum safely representable SNP name is 511 bytes, plus NUL.** This is a buffer-imposed limit, not a validated acceptance policy or a count of Unicode characters. Neither the R wrappers nor the native parser enforce it. Oversized tokens are copied rather than deliberately rejected; undefined behavior means there is no meaningful maximum reliably “accepted” length. A single read can contain up to 999 bytes before NUL; even a token filling that chunk is copied before missing-column checks run. Longer physical lines are split into chunks without explicit line-length validation.

`strcpy(ref, tmp)` at lines 124 and 134 also has no local check. Those copies are safe only if the preceding copy established a terminated name of at most 511 bytes. Once `tmp` has overflowed, subsequent comparisons/copies cannot be assumed safe.

### Earlier overflow in the counting pass

`nb_ind_ancestrymap` (`ancestrymap.c:202–233`) first reads into another 1000-byte line buffer, tokenizes it, and executes `strcpy(tmp, token)` at line 217, where `tmp` is again 512 bytes. Both public conversion functions call this before `read_ancestrymap`. Thus an oversized name on the first record can fail here before reaching line 160. An oversized name on a later locus can instead reach the reported parsing copy; the counting pass compares subsequent tokens without copying them into `tmp`.

The same function ignores the result of `fgets` before tokenizing, including at EOF, and does not check for a NULL token before `strcpy` or `strcmp`. Empty/blank input and EOF handling are additional malformed-input risks; the exact outcome depends on the buffer contents and where the first change of locus occurs. These must not be mistaken for deliberate length rejection.

## 2. R-to-native execution path

1. `NAMESPACE:8` loads the LEA shared library with `.registration=TRUE`; lines 31–32 export the ANCESTRYMAP conversion functions.
2. `R/ancestrymap2geno.R:2–30` checks argument type and the `.ancestrymap` filename extension, chooses an output filename when omitted, and optionally returns early when `force=FALSE` and both files exist. It never reads or validates SNP names.
3. Its `.C("R_ancestrymap2geno", ...)` call at line 24 passes input/output character values and two integer values for dimensions.
4. `src/R_init_LEA.c:21–24` registers `R_ancestrymap2geno` with four arguments (`STRSXP`, `STRSXP`, `INTSXP`, `INTSXP`); `R_init_LEA` calls `R_registerRoutines` at line 78. Registration does not validate file contents.
5. `src/R_ancestrymap2geno.c:29–34` dereferences the filename arguments and calls `ancestrymap2geno`.
6. `src/convert/ancestrymap.c:35–62`: `nb_ind_ancestrymap` → `nb_lines(input_file, 1000)` → dimension calculation/allocation → `read_ancestrymap` → `read_line_ancestrymap` → `write_geno`. The first unchecked copy is in the counting step; the reported copy is in the reading step.

`ancestrymap2lfmm` follows the equivalent R/native path and shares both vulnerable parsing functions (`ancestrymap.c:66–92`). `test_input_file` can invoke either conversion implicitly for higher-level analysis functions when conversion is needed (`R/test_input_file.R:29,39`; callers include `sNMF`, `pca`, `LFMM`, and `createDataSet`). The user need only import an affected file, not call `.C` manually. Existing output files can bypass conversion under `force=FALSE`; that is not input validation.

## 3. Other string operations and format coverage

| Location | Finding | Reachability / evidence |
| --- | --- | --- |
| `ancestrymap.c:160,217` | Unbounded SNP-name copies into 512-byte stack objects. | Both ANCESTRYMAP conversions; independently demonstrated with ASan. |
| `ancestrymap.c:124,134` | Unbounded copies between 512-byte name arrays. | Safe only after enforcing the name invariant at the initial copy sites. |
| `vcf2geno.c:119–120,262,269` | Nine metadata destinations each have 512 bytes; `read_cnv_info` uses unchecked `strcpy` for all nine VCF fixed columns. Any such field longer than 511 bytes can overflow its heap allocation, including ID, INFO, REF, or ALT. Filtering multi-base REF/ALT occurs after copying. | `R_vcf2geno` → `vcf2geno` → `read_cnv_info`. ASan independently confirms the first-field overflow. `vcf2lfmm` calls `vcf2geno` first. Long VCF fields can be legitimate format input. |
| `register_convert.c:37–42,50–52` | Copies command-line paths and appends suffixes without destination capacities. | Unsafe helper API for undersized caller buffers. Repository search found no calls to `analyse_param_convert`; do not claim a demonstrated R-path overflow through this helper. |
| `ped.c` | Uses `strtok`, checks token presence, and consumes allele characters directly; no name-token copy into a smaller fixed string buffer. | The reported SNP-name copy defect was not found in PED conversion. Its fixed line-size heuristic can split records with long IDs; this is a separate parsing concern. |
| `geno.c`, `geno2lfmm.c`, `lfmm2geno.c` | No analogous unbounded string copy in these files. GENO parses characters; LFMM conversion delegates numeric reading to `io_data_int.c`. | No evidence that a long SNP name triggers this particular defect in GENO/LFMM, which do not carry these name fields. This is not a general memory-safety certification. |

The string-operation review covered `strcpy`, `strcat`, tokenization, `strlen`, bounded line reads, and formatted I/O. No additional `sprintf`, `scanf`/`sscanf`, `gets`, or `memcpy` sites were found in `src/convert/`. Formatted output uses constant format strings; input strings used as `%s` arguments are not format-string injection sites.

### Related reachable wrapper risk outside `src/convert/`

`src/R_vcf2geno.c:34–43` constructs paths in stack buffers `snp_bp_file[512]` and `removed_bp_file[512]`, copying `remove_ext(*R_output_file, '.', '/')` and appending `.vcfsnp` / `.removed`. There is no capacity check. If the stem has length L, the writes require L+8 and L+9 bytes respectively, so both are safe only for L ≤ 503. This is a separate statically established risk, not dynamically tested here. Long nested output paths can reach it before opening the files. `remove_ext` itself allocates `strlen(input)+1` for its copy (`src/io/io_tools.c:44–77`); the smaller wrapper buffers are the problem.

### Additional issues observed, outside the minimal name fix

* `fill_line_vcf` writes `allele[i]` for every sample token without first checking `i < N`; `test_column` runs after the writes. Excess sample columns can therefore cause an out-of-bounds heap write. This is clear from source, but was not dynamically exercised in this audit.
* VCF header parsing assumes successful reads/non-NULL tokens; malformed files can reach invalid comparisons. ANCESTRYMAP has the analogous counting-pass read/token issues described above.
* Dimension arithmetic, unchecked allocations, and line-size heuristics deserve separate hardening. For example, `nb_lines(input_file, 1000)` reads chunks of up to 20009 bytes (`io_tools.c:138`), whereas ANCESTRYMAP parsing reads only 999 bytes at a time. Neither guarantees counting physical records consistently for oversized lines.

These findings should be tracked separately; an ANCESTRYMAP name check does not remedy VCF or establish that every conversion input is safe.

## 4. Dynamic verification performed

Environment: macOS arm64, Apple clang 21.0.0. The temporary harness and logs are under `/private/tmp/lea-ancestrymap-audit/` (not committed; temporary artifacts may later disappear).

The harness compiles the actual, unchanged `src/convert/ancestrymap.c` and `src/convert/vcf2geno.c`. A temporary `R.h` declares `Rprintf` and `error`; the harness implements logging/error exit and file opening. Unused conversion dependencies are stubs. It directly calls `read_line_ancestrymap`, `nb_ind_ancestrymap`, and `read_cnv_info` in separate processes. The line-parser destination is a real 512-byte stack array, matching the production caller; VCF destinations are nine separate 512-byte heap allocations, matching production. These are isolated native tests, not an end-to-end run inside R.

Build used, from repository root:

```sh
clang -g -O1 -fsanitize=address -fno-omit-frame-pointer \
  -I/private/tmp/lea-ancestrymap-audit -Isrc \
  /private/tmp/lea-ancestrymap-audit/harness.c \
  src/convert/ancestrymap.c src/convert/vcf2geno.c \
  -o /private/tmp/lea-ancestrymap-audit/harness
```

Inputs contain ASCII `A` repeated to the requested byte length. Line mode appends `" ind 1\n"`; count mode writes that record followed by `"B ind 0\n"`, ensuring a locus change without depending on EOF behavior. VCF mode appends `" 1 id A G . PASS . GT 0/1\n"` to the long first field. Run the retained harness, for example:

```sh
/private/tmp/lea-ancestrymap-audit/harness line 600
/private/tmp/lea-ancestrymap-audit/harness count 512 /private/tmp/lea-ancestrymap-audit/input.ancestrymap
/private/tmp/lea-ancestrymap-audit/harness vcf 600
```

| Routine / mode | 511 bytes | 512 bytes | 600 bytes |
| --- | --- | --- | --- |
| `read_line_ancestrymap` / `line` | Success, name preserved, allele 1 | ASan stack-buffer-overflow; write size 513 | ASan stack-buffer-overflow; write size 601 |
| `nb_ind_ancestrymap` / `count` | Success, N=1 | ASan stack-buffer-overflow; write size 513 | ASan `strcpy-param-overlap`, with destination range exceeding `tmp[512]` |
| `read_cnv_info` / `vcf` | Success | ASan heap-buffer-overflow; write size 513 | ASan heap-buffer-overflow; write size 601 |

All failing cases terminated via ASan (SIGABRT). The different counting-pass diagnostic at 600 bytes reflects this compiler's stack layout: the oversized destination range overlaps the source line buffer. It does not make that copy safe. The line-mode 600-byte result independently reproduces the reported 601-byte write.

## 5. Severity and consequences

Treat this as a high-priority native memory-safety defect (stack-based buffer overflow, CWE-121; the VCF metadata counterpart is heap-based, CWE-122). The demonstrated result is an out-of-bounds copy and sanitizer termination. Ordinary unsanitized builds have undefined behavior: crashes, corrupted local data, and unreliable conversion results are plausible consequences. The error affects the hosting R process and can terminate a research job or a service that imports uploaded data.

The attacker-controlled surface is file content processed by a conversion. No network listener is introduced by these routines, and an attacker must cause a user or service to process crafted input. In automated import services that may require no interactive action. Arbitrary code execution is a possible general consequence of native memory corruption, but no instruction-pointer control, exploit chain, information disclosure, privilege escalation, or code execution was demonstrated. Stack layout, compiler hardening, platform mitigations, and deployment context affect exploitability. A specific CVSS score or assertion of remote code execution would exceed the available evidence.

## 6. Proposed minimal correction — not implemented

Preserve all names of 1–511 bytes exactly and reject names of 512 bytes or more before copying. Do not use `strncpy`, silent truncation, or merely enlarge the buffers: truncation can merge distinct locus names, and enlargement leaves the unchecked contract unresolved.

A small coordinated correction should:

1. Define one shared internal capacity constant of 512 bytes for ANCESTRYMAP names and use it for `tmp` and `ref` declarations and checks. Make the destination contract explicit in the header. A small checked-copy helper can share the rule; alternatively pass the destination size into `read_line_ancestrymap` and update its declaration/caller.
2. In `read_line_ancestrymap`, after confirming the first token exists and **before** line 160, reject `strlen(token) >= capacity` with an R error identifying the input file, individual/locus position, and maximum of 511 bytes. Copy the full token and NUL only after validation (`memcpy(name, token, length + 1)` or checked `strcpy`).
3. Apply the same validation before line 217 in `nb_ind_ancestrymap`. Validate first-column tokens in its scan as well, so counting never treats an oversized token as valid. Check successful `fgets` and non-NULL tokens before using them. Handle EOF explicitly rather than tokenizing stale data; preserve the correct individual count for a valid final single-locus group.
4. Keep `ref` copies subject to the established 511-byte invariant. Both ANCESTRYMAP entry points must share the same policy; an R-only check is insufficient because native routines can be called directly.

The two copy-site length guards are the essential security change. The closely related read/token checks are needed to make rejection and regression testing dependable. Robust physical-line-length detection, cleanup across R errors, allocation/dimension checks, and the independently affected VCF/path routines should be explicit follow-up work rather than silently expanding this patch. Resource cleanup should close an opened input before raising an error where practical; a larger move to R-managed allocation/unwind cleanup is separate work.

## 7. Proposed regression tests — not added

Use temporary input/output files and isolated processes for cases that can crash the unfixed package. Add tests only after approval. Expected results apply to the corrected implementation.

| Case | Expected assertion |
| --- | --- |
| Name lengths 1, 510, 511 | Both ANCESTRYMAP converters succeed; compare exact expected GENO/LFMM output and dimensions. |
| Lengths 512, 513, 600 | Both converters raise an R error mentioning oversized SNP name and the 511-byte maximum; no output conversion is produced and no crash occurs. |
| Oversized name on first record | Exercises counting-pass guard, not just line parser. |
| Oversized name on second locus, after a short first locus | Exercises normal read path and line-parser guard; include a direct parser test to ensure this guard is covered even if counting rejects earlier. |
| Distinct 511-byte names sharing their first 510 bytes | Remain separate loci; output order/count is unchanged. |
| Two oversized names sharing their first 511 bytes | Reject; never merge names by truncation. |
| Spaces, tabs, LF, CRLF; optional final newline | Preserve existing valid parsing; limit is measured in token bytes. |
| Multibyte names with 510/511/512 encoded bytes | Enforce bytes including room for NUL, irrespective of displayed character count. |
| Empty file, blank first record, missing fields, EOF after one locus | Controlled errors for malformed input; correct count/output for a valid one-locus file, without NULL/stale-token access. |
| Tokens of 999 and 1000+ bytes; very long physical records | Controlled rejection, no overflow, no accidental acceptance of a split token. Track unrelated long-line limitations separately. |
| Existing short-name fixture with multiple individuals/loci and genotypes 0/1/2/9 | Exact output preservation for both conversions. |

A useful deterministic fixture has two individuals per locus: `name1 ind1 0`, `name1 ind2 1`, `name2 ind1 2`, `name2 ind2 9`. Expected GENO content is `01\n29\n`; expected LFMM rows contain `0 2` and `1 9`. Generate name strings programmatically and confirm their byte lengths. Check that an error does not overwrite a pre-existing output file. Do not rely on output from a crashed process as a test oracle.

For sanitizer verification:

* Rebuild the native parser harness with `-g -O1 -fsanitize=address -fno-omit-frame-pointer` and rerun every boundary case in a separate process. Verify correct success/error behavior **and** absence of ASan diagnostics. A nonzero exit alone is insufficient: the vulnerable code also exits nonzero when ASan aborts.
* Add explicit tests for both native copy sites, plus the `read_ancestrymap` caller's actual buffers. Retain pre-fix sanitizer evidence for the 600-byte reproduction; after correction the same input must reach the intended size error without memory corruption.
* In a separate temporary package copy and temporary R library, compile and link all package native code with ASan, using a compatible sanitizer-enabled R or the platform-appropriate sanitizer runtime setup. Run both exported R functions and the registered `.C` entry points. Use fresh builds to avoid reusing existing unsanitized `.o` files. Do not install over the user's normal package or rebuild in the working source tree.
* Keep VCF regressions separate: lengths 511/512/600 in each of the nine fields, oversized sample counts, and wrapper path-stem boundaries 503/504/505. Those require their own fixes; the proposed ANCESTRYMAP patch is not expected to make them pass.

Full-package R/ASan integration and proposed post-fix tests were not run in this audit. No fix has been implemented; approval is required before implementation begins.
