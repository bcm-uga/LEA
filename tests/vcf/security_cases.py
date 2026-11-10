"""Deterministic inputs and post-fix contracts, shared by isolated runners."""
HEADER = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
META = ('##fileformat=VCFv4.2\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##INFO=<ID=K,Number=1,Type=String,Description="Long annotation">\n')

def document(fields=None, samples=None, n=1, newline=True):
    fields = fields if fields is not None else ['1', '1', 'id', 'A', 'C', '.', 'PASS', '.', 'GT']
    samples = samples if samples is not None else ['0/1']
    text = META + '\t'.join(HEADER + ['S'+str(i+1) for i in range(n)]) + '\n'
    text += '\t'.join(fields + samples)
    return text + ('\n' if newline else '')

def success(name, fields=None, samples=None, n=1, genotype='1\n', **kwargs):
    fields = fields if fields is not None else ['1','1','id','A','C','.','PASS','.','GT']
    removed = len(fields[3]) > 1 or len(fields[4]) > 1
    return dict(name=name, text=document(fields, samples, n), expectation='ok', n=n,
                m=0 if removed else len(genotype.splitlines()),
                geno='' if removed else genotype,
                vcfsnp='' if removed else ' '.join(fields)+' \n',
                removed=' '.join(fields)+' REMOVED\n' if removed else '', **kwargs)

def cases():
    items = []
    for index, label in enumerate(('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT')):
        for size in (511,512,600):
            fields = ['1','1','id','A','C','.','PASS','.','GT']
            fields[index] = ('K='+'a'*(size-2) if label == 'INFO' else
                             'GT:'+'X'*(size-3) if label == 'FORMAT' else
                             'A'*size if label in ('REF','ALT') else
                             '1'*size if label in ('POS','QUAL') else 'x'*size)
            c = success(f'metadata-{label}-{size}', fields,
                        ['0/1:.'] if label == 'FORMAT' else ['0/1'])
            # Numeric stress tokens need not be accepted if validation tightens;
            # they must never corrupt memory. All semantically valid long fields
            # require complete successful output, including removed REF/ALT.
            if label in ('POS','QUAL'):
                c['expectation'] = 'safe'
            if label in ('REF','ALT'):
                c['r_modes'] = ['geno']  # Empty-LFMM handling is characterized separately.
            items.append(c)
    fields = ['1','1','id','A','C','.','PASS','K='+'a'*5000,'GT']
    items.append(success('long-info-5002', fields))
    fields = ['1','1','id','A','C','.','PASS','.','GT:XX']
    items.append(success('long-physical-sample-2004', fields, ['0/1:'+'a'*2000]))
    items.append(success('no-final-newline'))
    items[-1]['text'] = items[-1]['text'].rstrip('\n')
    items.append(success('header-long-name'))
    items[-1]['text'] = items[-1]['text'].replace('\tS1\n','\t'+'S'*5000+'\n')
    for name,samples,n in [('fewer-samples',['0/1'],2),('extra-samples',['0/1','0/0'],1)]:
        items.append(dict(name=name,text=document(samples=samples,n=n),expectation='error'))
        fields=['1','1','id','AT','A','.','PASS','.','GT']
        items.append(dict(name='removed-'+name,text=document(fields,samples,n),expectation='error'))
    for value in ('0/', '/1', '0x1', '0/10', '0/1/0', '2', '-1', ''):
        items.append(dict(name='bad-gt-'+(value.replace('/','_') or 'empty'),
                          text=document(samples=[value]),expectation='error'))
    for fmt,gt,label in [('DP','1','no-GT'),('DP:GT','5:0/1','GT-not-first'),('GT:DP','0:','empty-subfield')]:
        fields=['1','1','id','A','C','.','PASS','.',fmt]
        items.append(dict(name='format-'+label,text=document(fields,[gt]),expectation='error'))
    # Empty ID is currently collapsed and causes later columns to shift.
    fields=['1','1','','A','C','.','PASS','.','GT']
    items.append(dict(name='empty-tab-ID',text=document(fields),expectation='error'))
    # An empty sample between nonempty samples must not be collapsed.
    items.append(dict(name='empty-tab-sample',text=document(samples=['0/1','','0/0'],n=3),expectation='error'))
    items.append(dict(name='empty-tab-sample-compensated',
                      text=document(samples=['0/1','','0/0'],n=2),expectation='error'))
    items.append(dict(name='empty-tab-ID-compensated',
                      text=document(fields,['0/1','0/0']),expectation='error'))
    items.append(dict(name='missing-fixed-columns',text=document(['1','1','id']),expectation='error'))
    items.append(dict(name='empty-file',text='',expectation='error'))
    items.append(dict(name='metadata-only',text=META,expectation='error'))
    items.append(dict(name='header-only',text=META+'\t'.join(HEADER+['S1'])+'\n',
                      expectation='ok',n=1,m=0,geno='',vcfsnp='',removed='',r_modes=['geno']))
    # Separate observation: existing geno2lfmm cannot safely convert zero loci.
    items.append(dict(name='header-only-lfmm',text=items[-1]['text'],expectation='safe', n=1, m=0, geno='',vcfsnp='',removed='',
                      r_modes=['lfmm'],native=False))
    for size in (503,504,505,511,512,600):
        items.append(success('path-stem-'+str(size),path_length=size,r_modes=['geno']))
    return items


def verdict(case, returncode, stdout, stderr, timed_out, files, native=False):
    """Require a success contract or a real controlled error, never just nonzero."""
    log = stdout + stderr
    if timed_out:
        return False, 'timeout'
    if 'AddressSanitizer' in log:
        kind = next((s.strip() for s in log.splitlines() if 'ERROR: AddressSanitizer' in s), 'ASan diagnostic')
        return False, 'ASan-'+kind.split('ERROR: AddressSanitizer:',1)[-1].strip().split()[0]
    controlled = returncode == 2 and ('CONTROLLED_ERROR:' if native else 'R_ERROR:') in log
    if returncode != 0:
        if controlled and case['expectation'] in ('error','safe'):
            return True, 'controlled-error'
        return False, ('legacy-null-error' if 'LEGACY_NULL_ERROR' in log else
                       f'signal-{-returncode}' if returncode < 0 else f'exit-{returncode}')
    if case['expectation'] == 'error':
        return False, 'unexpected-acceptance'
    for ext in ('geno','vcfsnp','removed'):
        if files.get(ext) != case[ext]:
            return False, 'wrong-'+ext
    if native and f'VCF_COUNTS:{case["n"]}:{case["m"]}' not in stdout:
        return False, 'wrong-counts'
    return True, 'success'
