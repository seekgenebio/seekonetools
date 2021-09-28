import os
import json
from functools import partial
from collections import defaultdict, Counter
import dnaio
from cutadapt import adapters
from ..utils.pipeline import Pipeline
from .helper import prepare_funcs
from ..utils._version import __version__


class Stat:
    def __init__(self):
        self.data = {}

    def update(self, **d):
        if not self.data:
            self.data = d
        else:
            for k, v in d.items():
                self.data[k] += v

    @staticmethod
    def _sort_gc(d):
        idx_max = max([k[0] for k in d])
        return {
            b: [d.get((i, b), 0) for i in range(idx_max)] for b in 'ATCGN'
        }

    @staticmethod
    def _sort_q(d, phred=33):
        idx_max = max([k[0] for k in d])
        q_max = max([ord(k[1])-phred for k in d])
        return {
            i: [d.get((i, chr(q+phred)), 0) for q in range(q_max+1)] for i in range(idx_max)
        }      

    def save(self, path='summary.json'):
        tmp = {'__version__': __version__}
        for k in self.data:
            if k.endswith('_gc'):
                tmp[k] = self._sort_gc(self.data[k])
            elif k.endswith('_q'):
                tmp[k] = self._sort_q(self.data[k])
            else:
                tmp[k] = dict(self.data[k])
        with open(path, 'w') as fh:
            json.dump(tmp, fh, indent=4)

def process_barcode(fq1, fq2, fq_out, shift, shift_pattern,
                r1_structure, funcs,  minlen=50):
    Barcode_GC_Counter = Counter()
    UMI_GC_Counter = Counter()
    R2_GC_Counter = Counter()
    Barcode_Q_Counter = Counter()
    UMI_Q_Counter = Counter()
    R2_Q_Counter = Counter()
    stat_Dict = defaultdict(int)

    start_pos = 0
    end_pos = 0
    flag = 0
    shift_pos = 0
    polyA_filter = adapters.BackAdapter(sequence='AAAAAAAAAAAAAAA')

    outfh = dnaio.open(fq_out, fileformat='fastq', mode='w')
    fh = dnaio.open(file1=fq1, file2=fq2, fileformat='fastq', mode='r')
    for r1, r2 in fh:
        stat_Dict['total'] += 1
        sequence = r1.sequence
        qualities = r1.qualities
        start_pos = 0
        end_pos = 0
        flag = 0
        old_seqs = defaultdict(list)
        new_seqs = defaultdict(list)
        seq_quals = defaultdict(str)
        if shift:
            shift_pos = sequence[:7].find(shift_pattern)
            if shift_pos < 0:
                stat_Dict['no_anchor'] += 1
                continue
            else:
                start_pos = shift_pos + 1

        for (code, n), func in zip(r1_structure, funcs):
            end_pos = start_pos + n
            seq = sequence[start_pos:end_pos]
            quals = qualities[start_pos:end_pos]
            
            old_seqs[code].append(seq)
            if func:
                # 处理需要做纠错的code
                base3 = sequence[end_pos:end_pos + 3]
                correct_code, indel, seqx = func(seq, seqafter=base3)

                if correct_code < 1:
                    # 0: 纠错无效，1: 无需纠错，2: 缺失纠错 3: 错配纠错，-1: 多种纠错可能
                    flag = 1
                    break
                else:
                    new_seqs[code].append(seqx)
                    seq_quals[code] += quals
                    start_pos = start_pos + n + indel
            else:
                # 处理不需要做纠错的code
                new_seqs[code].append(seq)
                seq_quals[code] += quals
                start_pos = start_pos + n

        if flag:
            stat_Dict[code] += 1
            continue

        stat_Dict['valid'] += 1
        barcode_base = ''.join(new_seqs['B'])
        umi_base = new_seqs['U'][0]

        # filter polyA
        match = polyA_filter.match_to(r2.sequence)
        if match:
            r2 = match.trimmed(r2)
            _len = len(r2)
            if _len < minlen:
                stat_Dict['too_short'] += 1
                continue
            else:
                stat_Dict['trimmed'] += 1

        # barcode_umi_:oldpart:_readID
        __barcode = ':'.join( ['' if n==o else o for n, o in zip(new_seqs['B'], old_seqs['B'])] )
        r2.name = '_'.join([barcode_base, umi_base, __barcode, r2.name])
        outfh.write(r2)

        Barcode_GC_Counter.update(enumerate(barcode_base))
        UMI_GC_Counter.update(enumerate(umi_base))
        R2_GC_Counter.update(enumerate(r2.sequence))
        Barcode_Q_Counter.update(enumerate(seq_quals['B']))
        UMI_Q_Counter.update(enumerate(seq_quals['U']))
        R2_Q_Counter.update(enumerate(r2.qualities))

    outfh.close()
    fh.close()
    return {
            'stat': Counter(stat_Dict),
            'barcode_gc':Barcode_GC_Counter,
            'umi_gc': UMI_GC_Counter,
            'r2_gc': R2_GC_Counter,
            'barcode_q': Barcode_Q_Counter,
            'umi_q': UMI_Q_Counter,
            'r2_q': R2_Q_Counter
        }


def barcode(fq1:list, fq2:list, samplename: str, outdir:str,
            barcode:list=[], shift:str=True, shift_pattern:str='A',
            structure:str='B8L15B8L15B8U12T15', linker: list=[],
            B:tuple=(1,0), L:tuple=(1,0), core:int=4, logger:object=None,
            **kwargs):

    logger.info('cellbarcode started!')
    outdir1 = os.path.join(outdir, 'step1')
    os.makedirs(outdir1, exist_ok=True)

    r1_structure, funcs = prepare_funcs(
                            structure=structure,
                            _barcode=barcode,
                            _linker=linker,
                            B=B, L=L)
    logger.info('prepare done!')
    # pre precess? check structure?

    # start worker processes
    worker_func = partial(
                    process_barcode,
                    shift=shift,
                    shift_pattern=shift_pattern,
                    r1_structure=r1_structure,
                    funcs=funcs
                    )

    stat = Stat()
    pipeline = Pipeline(
        func = worker_func,
        fq1 = fq1,
        fq2 = fq2,
        fqout1=os.path.join(outdir1, f'{samplename}_2.fq.gz'),
        stat=stat,
        core=core,
    )
    pipeline.run()
    pipeline.stat.save(os.path.join(outdir, f'{samplename}_summary.json'))
