#!/usr/bin/env python3
"""Instrumented early-control-byte check under ASan/UBSan; build only in /tmp."""
import os
from pathlib import Path
import subprocess
import tempfile

repo = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="lea-vcf-control-") as temporary:
    root = Path(temporary)
    (root / "R.h").write_text("void Rprintf(const char *, ...);\n_Noreturn void error(const char *, ...);\n")
    binary = root / "control-bytes"
    subprocess.run([os.environ.get("CC", "clang"), "-g", "-O1", "-std=c11",
                    "-fsanitize=address,undefined", "-fno-sanitize-recover=undefined",
                    "-fno-omit-frame-pointer", "-I" + str(root), "-I" + str(repo / "src"),
                    str(repo / "tests/native/vcf-control-bytes.c"), "-o", str(binary)], check=True)
    env = dict(os.environ, ASAN_OPTIONS="detect_leaks=0:halt_on_error=1")
    result = subprocess.run([str(binary)], capture_output=True, text=True, env=env, cwd=root, timeout=10)
    log = result.stdout + result.stderr
    assert result.returncode == 0 and "AddressSanitizer" not in log and "runtime error:" not in log, log
    print(log, end="")
