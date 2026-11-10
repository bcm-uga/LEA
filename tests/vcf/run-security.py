#!/usr/bin/env python3
"""Isolated VCF security contracts (intentionally failing before the fix).
Default: enforce post-fix contracts and exit 1 on any failure.
--characterize: record current behavior, report failures, exit 0 (NOT a pass).
--native: compile/run AddressSanitizer driver; otherwise use Rscript children.
--library /path: required R library for real R integration tests.
--report /outside/repo/results.json: optional full observations and child logs.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import sys

# Keep import caches out of the repository.
sys.dont_write_bytecode = True

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
spec = importlib.util.spec_from_file_location('vcf_cases', HERE/'security_cases.py')
cases_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cases_module)

def long_output(root, size):
    # Length of complete absolute stem, not merely its final component.
    parent = root
    while size - len(str(parent)) - 1 > 80:
        parent = parent / ('d' * min(60, size-len(str(parent))-1-40))
    parent.mkdir(parents=True, exist_ok=True)
    remaining = size - len(str(parent)) - 1
    if remaining < 1:
        raise RuntimeError('Temporary root too long for requested stem boundary')
    stem = parent / ('x'*remaining)
    assert len(str(stem)) == size
    return Path(str(stem)+'.geno')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', action='store_true')
    parser.add_argument('--library')
    parser.add_argument('--characterize', action='store_true')
    parser.add_argument('--report', type=Path)
    parser.add_argument('--timeout', type=float, default=3)
    args = parser.parse_args()
    if not args.native and not args.library:
        parser.error('--library is required for the R runner')
    observations = []
    with tempfile.TemporaryDirectory(prefix='lea-vcf-security-') as tmp:
        root = Path(tmp)
        env = dict(os.environ)
        if args.native:
            (root/'R.h').write_text('void Rprintf(const char *, ...);\n_Noreturn void error(const char *, ...);\n')
            sources = ['tests/native/vcf-security.c','src/R_vcf2geno.c','src/convert/vcf2geno.c',
                       'src/convert/geno.c','src/convert/register_convert.c','src/io/io_tools.c',
                       'src/io/io_error.c','src/matrix/error_matrix.c']
            binary = root/'vcf-asan'
            subprocess.run([os.environ.get('CC','clang'),'-g','-O1','-std=c11','-fsanitize=address',
                            '-fno-omit-frame-pointer','-I'+str(root),'-I'+str(REPO/'src'),
                            *[str(REPO/p) for p in sources],'-o',str(binary)],check=True)
            env['ASAN_OPTIONS'] = 'detect_leaks=0:halt_on_error=1'
        else:
            env['R_LIBS'] = str(Path(args.library).resolve())
        for case in cases_module.cases():
            if args.native and not case.get('native', True):
                continue
            modes = ['native'] if args.native else case.get('r_modes',['geno','lfmm'])
            for mode in modes:
                work = root/(case['name']+'-'+mode); work.mkdir()
                input_path = work/'input.vcf'; input_path.write_bytes(case['text'].encode())
                output = long_output(work,case['path_length']) if 'path_length' in case else work/'output.geno'
                cmd = ([str(binary),str(input_path),str(output)] if args.native else
                       ['Rscript',str(HERE/'isolated-security.R'),mode,str(input_path),str(output)])
                timed_out = False
                try:
                    result = subprocess.run(cmd,cwd=work,env=env,capture_output=True,text=True,timeout=args.timeout)
                    status,stdout,stderr = result.returncode,result.stdout,result.stderr
                except subprocess.TimeoutExpired as e:
                    timed_out = True; status = None
                    stdout = e.stdout or b''; stderr = e.stderr or b''
                    if isinstance(stdout,bytes): stdout=stdout.decode(errors='replace')
                    if isinstance(stderr,bytes): stderr=stderr.decode(errors='replace')
                # LFMM wrapper uses the input-derived intermediate GENO path.
                stem = input_path.with_suffix('') if mode == 'lfmm' else output.with_suffix('')
                files = {ext: Path(str(stem)+'.'+ext).read_text() if Path(str(stem)+'.'+ext).exists() else None
                         for ext in ('geno','vcfsnp','removed','lfmm')}
                passed,reason = cases_module.verdict(case,status,stdout,stderr,timed_out,files,args.native)
                if passed and mode == 'lfmm' and case['expectation'] != 'error' and status == 0:
                    # Expected matrix orientation: individuals in rows, loci in columns.
                    rows=case['geno'].splitlines()
                    expected=(''.join(' '.join(row[i] for row in rows)+'\n' for i in range(case['n']))
                              if rows else '')
                    if files['lfmm'] != expected: passed,reason=False,'wrong-lfmm'
                observations.append(dict(case=case['name'],mode=mode,expectation=case['expectation'],
                                         passed=passed,reason=reason,returncode=status,
                                         stdout=stdout,stderr=stderr,files=files))
                print(('PASS' if passed else 'FAIL')+' '+case['name']+'/'+mode+': '+reason)
    if args.report:
        args.report.write_text(json.dumps(observations,indent=2)+'\n')
    failures = sum(not item['passed'] for item in observations)
    print(f'{len(observations)} isolated cases; {failures} post-fix contract failures')
    if args.characterize:
        print('CHARACTERIZATION ONLY: failures are observations, not waived security tests.')
    return 0 if args.characterize or failures == 0 else 1

if __name__ == '__main__':
    raise SystemExit(main())
