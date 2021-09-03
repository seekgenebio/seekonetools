import re
import sys
from subprocess import run
from collections import defaultdict
from xopen import xopen


def gtfstat(gtf, feature="gene", key='gene_biotype'):
    '''
    gtf summary
    '''
    with xopen(gtf) as fh:
        summary = defaultdict(int)
        for line in fh:
            if line.startswith('#'): continue
            if not line.strip(): continue
            tmp = line.strip().split('\t')
            if tmp[2] == feature:
                tmp_list = [s.strip() for s in tmp[-1].strip().strip(';').split(';')]
                tmp_dict = dict([l.replace('"', '').split(' ', 1) for l in tmp_list])
                summary[tmp_dict[key]] += 1
    for k, v in sorted(summary.items(), key = lambda x: x[1]):
        sys.stdout.write(f'{k}\t{v}\n')

def gtfparse(gtf):
    '''
    parse gtf
    '''
    with xopen(gtf) as fh:    
        record = []
        for l in fh:
            if not l.strip():
                continue
            if l.startswith('#'):
                record.append([l.strip()])
                continue
            tmp = l.strip().split('\t')
            if tmp[2] == 'gene':
                yield record
                record = []
            else:
                record.append(tmp)
        yield record
        
def gtffilter(gtf, biotype, key='gene_biotype', filtered_gtf=None):
    '''
    filter gtf
    '''
    print(gtf)
    print(biotype)
    print(key)
    print(filtered_gtf)
    if not filtered_gtf:
        filtered_gtf = gtf.replace('.gtf', '.filtered.gtf')
    gtf_g = gtfparse(gtf)
    with xopen(filtered_gtf, 'w') as fhout:
        for r in gtf_g:
            gene_id = None
            gene_name = None
            for tmp in r:
                if tmp[0].startswith('#'):
                    fhout.write(f'{tmp[0]}\n')
                elif tmp[2]=='gene':
                    tmp_list = [s.strip() for s in tmp[-1].strip().strip(';').split(';')]
                    tmp_dict = dict([l.replace('"', '').split(' ') for l in tmp_list])
                    if tmp_dict[key] in biotype:
                        print("yes")
                        gene_id = tmp_dict['gene_id']
                        gene_name = tmp_dict.get('gene_name', gene_id)
                        if 'gene_name' not in tmp[-1]:
                            tmp[-1] = f'{tmp[-1].strip(";")}; gene_name "{gene_name}"'

                        if tmp[0]=='MT': 
                            if not re.match('MT', gene_name, re.I):
                                gene_name = f'MT-{gene_name}'
                                tmp[-1] = re.sub(r'gene_name ".*?"', 'gene_name "{gene_name}"')
                        fhout.write('{}\n'.format("\t".join(tmp)))
                    else:
                        break
                else:
                    if 'gene_name' not in tmp[-1]:
                        tmp[-1] = f'{tmp[-1].strip(";")}; gene_name "{gene_name}"'
                    fhout.write('{}\n'.format("\t".join(tmp)))


def mkref(fa, gtf, geomeDir, star_path='STAR'):
    cmd = (f'{star_path} --runMode genomeGenerate --runThreadN 16 '
           f'--genomeDir {geomeDir} --genomeFastaFiles {fa} '
           f'--genomeSAindexNbases 14   --genomeChrBinNbits 18 '
           f'--genomeSAsparseD 3   --limitGenomeGenerateRAM 17179869184 '
           f'--sjdbGTFfile {gtf}'
           )    
    run(cmd)


