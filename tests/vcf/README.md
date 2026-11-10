# VCF compatibility baseline and isolated security tests

Baseline: LEA 3.25.1, commit `c56f629`, before any VCF implementation changes.
Production files are unchanged. Build/install the baseline in a temporary source
copy and library; do not reuse existing compiled objects from the working tree.

## Safe default package test

```sh
R_LIBS=/path/to/baseline/library Rscript tests/vcf-compatibility.R
```

This is the only new top-level `.R` test. It runs safe fixtures in-process:
8 fixtures × 5 operations (two GENO wrappers, two LFMM wrappers, registered native
GENO call), plus 8 wrapper path/force scenarios. It compares exact output bytes,
returned path strings, native N/M counts, and captured stdout including removed
counts. Outputs/build artifacts go to temporary directories.

`fixtures/` contains:

- `bundled.vcf`: bytes from `data/example_vcf.txt.gz` with `#` prepended to its
  descriptive `CHROM` heading, producing the `#CHROM` heading used by LEA's
  documented conversion example. No variant records are changed.
- `diploid-tabs.vcf`: unphased/phased, completely/partially missing, retained and
  removed variants. `space-legacy.vcf`, `mixed-legacy.vcf`, `crlf.vcf`, and
  `no-final-newline.vcf` carry the same data with different separators/endings.
- `haploid-tabs.vcf`: haploid calls, missing calls, and FORMAT suffixes.
- `suffix-tabs.vcf`: diploid FORMAT suffixes and valid omitted trailing values.

`expected/<fixture>/output.{geno,lfmm,vcfsnp,removed}` are captured baseline bytes.
Sidecar trailing spaces and empty files are deliberate. Do not normalize them.
`expected/behavior.R` stores relative returned paths (including `./`), stdout,
native counts, and exact force-scenario file snapshots. These expectations do not
assert dangerous malformed-input behavior.

`baseline-manifest.json` records fixture/golden SHA-256 values and source/version
provenance. `baseline-source-md5.R` protects the explicit capture command against
running after implementation changes in the five relevant C/R source files.
To reproduce the capture **only against the unmodified baseline library**, run:

```sh
R_LIBS=/path/to/baseline/library Rscript tests/vcf/capture-baseline.R
```

Capture overwrites golden files by explicit request; normal tests never do so.
Never regenerate goldens from the security-fixed implementation merely to make a
regression pass. Compare with the recorded baseline and review any discrepancy.

## Security tests: isolated, intentionally red before the fix

```sh
python3 tests/native/run-vcf-asan.py
python3 tests/vcf/run-security.py --library /path/to/baseline/library
```

Both commands enforce the proposed safe behavior and exit **1** on current defects.
They are deliberately not invoked by the default R package test. Each case runs
in a separate process with a timeout. No dangerous input is converted in the main
R regression process. The R child is `isolated-security.R`; do not source it in
an interactive session. The native driver calls the actual registered wrapper,
links production code unchanged, and substitutes only R logging/error routines.

For observation-only runs:

```sh
python3 tests/native/run-vcf-asan.py --characterize --report /tmp/vcf-native.json
python3 tests/vcf/run-security.py --library /path/to/library --characterize --report /tmp/vcf-r.json
```

`--characterize` exits 0 after collecting observations even when contracts fail;
it is **not** a security pass. `--report` should point outside the repository.
`--timeout` sets the per-child deadline (default 3 seconds). Compiler/build
failures remain fatal. Python import caches are disabled for local case modules.
All compiler outputs and generated VCFs/paths are temporary.

Contracts in `security_cases.py`:

- `ok`: successful conversion, exact GENO/sidecar bytes, native counts, and LFMM
  orientation when exercised. Valid long metadata must succeed without truncation.
- `error`: controlled informative error; neither crashes, hangs, silent acceptance,
  nor a native `error(NULL)` stub exit counts as success.
- `safe`: either exact successful output or a controlled error. Used only for
  synthetic oversized numeric tokens and the unresolved header-only LFMM case.
  These do not require preservation of semantically invalid input.

The native suite checks address accesses with ASan. Leak checks are disabled
because the native R-error replacement exits; resource ownership is not validated
by that stub. Real R children separately exercise actual R error propagation.
Both runners inspect output content; a zero exit without correct files is not a
pass, and an arbitrary nonzero exit is not a controlled error.

`observed-characterization.tsv` records the current run's 60 native and 108 R
observations. It is documentation, **not** a list of accepted expected failures.
Native diagnoses are stable evidence; unsanitized corruption may manifest as
wrong output or different signals on another run/platform. Detailed traces stay
in external JSON/log files. See `VCF_SECURITY_BASELINE.md` for decisions and limits.
