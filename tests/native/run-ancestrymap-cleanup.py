#!/usr/bin/env python3
"""Fault-inject the production ANCESTRYMAP context under ASan/UBSan."""
import os
from pathlib import Path
import subprocess
import tempfile

repo = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="lea-ancestrymap-faults-") as tmp:
    root = Path(tmp)
    (root / "R.h").write_text("void Rprintf(const char *, ...);\n_Noreturn void error(const char *, ...);\n")
    fixture = root / "input.ancestrymap"
    fixture.write_text("a " + "I"*5000 + " 0 note\na j 3\nb i 9\nb j -1")
    binary = root / "faults"
    subprocess.run([os.environ.get("CC", "clang"), "-std=c11", "-g", "-O1",
                    "-fsanitize=address,undefined", "-fno-sanitize-recover=undefined",
                    "-fno-omit-frame-pointer", "-I"+str(root), "-I"+str(repo/"src"),
                    str(repo/"tests/native/ancestrymap-cleanup.c"), "-o", str(binary)], check=True)
    env = dict(os.environ, ASAN_OPTIONS="detect_leaks=0:halt_on_error=1")
    result = subprocess.run([str(binary), str(fixture), str(root/"output")], env=env, cwd=root,
                            capture_output=True, text=True, timeout=30)
    log = result.stdout + result.stderr
    assert result.returncode == 0 and "AddressSanitizer" not in log and "runtime error:" not in log, log
    print(log, end="")
