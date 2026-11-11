#!/usr/bin/env python3
"""Thin entry point for the isolated public-wrapper ASan security suite."""
from pathlib import Path
import subprocess
import sys
runner = Path(__file__).resolve().parents[1]/'vcf'/'run-security.py'
raise SystemExit(subprocess.call([sys.executable, str(runner), '--native', *sys.argv[1:]]))
