#!/usr/bin/env python3
"""Run from any directory: python3 tests/native/run-ancestrymap-asan.py.
Requires clang (or CC) and an AddressSanitizer runtime; builds only in /tmp.
"""
import os
from pathlib import Path
import subprocess
import tempfile

repo = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix="lea-ancestrymap-asan-") as directory:
    root = Path(directory)
    # Standalone replacements for the two R functions used by these sources.
    (root / "R.h").write_text(
        "void Rprintf(const char *, ...);\n"
        "_Noreturn void error(const char *, ...);\n")
    sources = ["tests/native/ancestrymap-security.c", "src/convert/ancestrymap.c",
               "src/convert/geno.c", "src/io/io_tools.c", "src/io/io_error.c",
               "src/io/io_data_int.c", "src/matrix/error_matrix.c"]
    executable = root / "ancestrymap-asan"
    subprocess.run([os.environ.get("CC", "clang"), "-g", "-O1", "-std=c11",
                    "-fsanitize=address", "-fno-omit-frame-pointer",
                    "-I" + str(root), "-I" + str(repo / "src"),
                    *[str(repo / s) for s in sources], "-lm", "-o", str(executable)],
                   check=True)
    checks = 0
    def execute(mode, argument, expected=None, message=None):
        global checks
        output = root / "output"
        output.write_text("sentinel\n")
        env = dict(os.environ)
        # R errors use nonlocal exits in production. This driver exits instead;
        # disable leak reporting, but keep all ASan memory-access checks active.
        env["ASAN_OPTIONS"] = "detect_leaks=0:halt_on_error=1"
        result = subprocess.run([str(executable), mode, str(argument), str(output)],
                                text=True, capture_output=True, env=env)
        log = result.stdout + result.stderr
        assert "AddressSanitizer" not in log, log
        assert result.returncode == (2 if message else 0), log
        if message:
            assert message in log, log
            assert output.read_text() == "sentinel\n"
        elif expected is not None:
            assert output.read_text() == expected, output.read_text()
        checks += 1

    fixture = root / "input.ancestrymap"
    for length in (1, 510, 511, 512, 513, 600):
        name = "a" * length
        message = "maximum is 511 bytes" if length >= 512 else None
        execute("line", length, message=message)
        for first in (True, False):
            a, b = (name, "b") if first else ("b", name)
            fixture.write_text(f"{a} i 0\n{a} j 1\n{b} i 2\n{b} j 9\n")
            for mode, expected in (("geno", "01\n29\n"), ("lfmm", "0 2\n1 9\n")):
                execute(mode, fixture, expected, message)
        fixture.write_text(f"b i 0\nb j 1\nc i 2\nc j 9\n{name} i 0\n{name} j 1\n")
        execute("geno", fixture, "01\n29\n01\n", message)
        execute("lfmm", fixture, "0 2 0\n1 9 1\n", message)
    for text, geno, lfmm, message in (
        ("a i 0\na j 1\n", "01\n", "0\n1\n", None),
        ("a i 0\na j 1", "01\n", "0\n1\n", None),
        ("a i 2", "2\n", "2\n", None),
        ("", None, None, "Empty ANCESTRYMAP"),
        ("\n", None, None, "missing SNP name"),
        ("a", None, None, "missing individual identifier"),
        ("a i", None, None, "missing genotype"),
        ("a i 0\nb i", None, None, "missing genotype"),
        ("a i 0\na j 1\nb i 2", None, None, "incorrect number of lines"),
        ("a i 0\na j 1\nb i 2\nc i 1\nc j 2\nc k 9", None, None,
         "individuals instead of"),
        ("a i 0\n\n", None, None, "missing SNP name"),
        ("a i 0\na j 1\nb i 2\nb j 9\nb k 0\nc i 1\n", None, None,
         "Invalid ANCESTRYMAP group sizes"),
    ):
        fixture.write_text(text)
        execute("geno", fixture, geno, message)
        execute("lfmm", fixture, lfmm, message)
    print(f"ASan ANCESTRYMAP checks passed: {checks}; no sanitizer diagnostics")
