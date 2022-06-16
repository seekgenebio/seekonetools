import os
import json
from collections import defaultdict
from subprocess import run


def run_STAR(fq, genomeDir, prefix, core=4, star_path='STAR'):
    args = [
        star_path, '--runThreadN', core, '--limitOutSJcollapsed', 5000000,
        '--genomeDir', genomeDir, '--readFilesIn', fq, '--readFilesCommand',
        'zcat', '--outFileNamePrefix', prefix, '--outSAMtype', 'BAM',
        'Unsorted'
    ]
    args = [str(_) for _ in args]
    # logger.info(' '.join(args))
    call_info = run(args, check=False)
    return f'{prefix}Aligned.out.bam', f'{prefix}Log.final.out'

def run_bamsort(inbam, outbam, byname=False, clean=True, core=4, samtools_path='samtools'):
    args = [
        samtools_path, 'sort', '-O', 'BAM', '-@',  core, '-o', outbam, inbam
    ]
    if byname:
        args.insert(2, '-n')
    args = [str(_) for _ in args]
    call_info = run(args, check=True)
    if call_info.returncode == 0 & clean:
        os.remove(inbam)
    return outbam

def run_qualimap(bam, gtf, outdir, SC5P, s="f"):
    strand = { 
        'f': 'strand-specific-forward',
        'r': 'strand-specific-reverse',
        'non': 'non-strand-specific'
    }   
    if SC5P:
        s="r"
    # gtf.gz?
    if gtf.endswith('.gz'):
        gtf_file = os.path.join(outdir, 'tmp.gtf')
        args = [
            'gzip', '-dc', gtf, '>', gtf_file
        ]
        assert not os.system(' '.join(args)), f"failed to unzip {gtf}!!!"
    else:
        gtf_file = gtf
    args = [
        'qualimap', 'rnaseq', '-outformat', 'PDF', '-outdir', outdir, '-bam',
        bam, '-gtf', gtf, '-p', strand[s], '--java-mem-size=8G'
    ]
    my_env = os.environ.copy()
    if 'DISPLAY' in my_env:
        del my_env['DISPLAY']
    run(args, check=False, env=my_env)
    return os.path.join(outdir, 'rnaseq_qc_results.txt')

def mapping_summary(STARLog, RnaSeqMetrics):
    summary = defaultdict()
    with open(STARLog, 'r') as fh:
        for line in fh:
            if 'Number of input reads' in line:
                summary['Number of input reads'] = int(
                    line.strip().split('\t')[-1])
            if 'Uniquely mapped reads number' in line:
                summary['Uniquely mapped reads number'] = int(
                    line.strip().split('\t')[-1])
            if 'Number of reads mapped to multiple loci' in line:
                summary['Number of reads mapped to multiple loci'] = int(
                    line.strip().split('\t')[-1])
    with open(RnaSeqMetrics, 'r') as fh:
        while True:
            line = fh.readline().strip()
            if line.startswith('total alignments'):
                summary['total alignments'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('reads aligned'):
                summary['reads aligned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('aligned to genes'):
                summary['aligned to genes'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('no feature assigned'):
                summary['no feature assigned'] = int(line.split()[-1].replace(
                    ',', ''))
            if line.startswith('exonic'):
                summary['exonic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intronic'):
                summary['intronic'] = int(line.split()[-2].replace(',', ''))
            if line.startswith('intergenic'):
                summary['intergenic'] = int(line.split()[-2].replace(',', ''))
                break
    return summary


def run_featureCounts(bam, gtf, samplename, outdir, region, SC5P, core=4, **kwargs):
    outcounts = os.path.join(outdir, 'counts.txt')
    s="1"
    if SC5P:
        s="2"
    args = [
        'featureCounts', '-T', core, '-t', region ,'-s', s, '-M', '-O', '-g', 'gene_id',
        '--fracOverlap', 0.5, '-a', gtf, '-o', outcounts, '-R', 'BAM', bam
    ]
    args = [str(_) for _ in args]
    run(args, check=True)
    featureCounts_bam = os.path.join(
        outdir, f'{samplename}_SortedByCoordinate.bam.featureCounts.bam')
    return featureCounts_bam

def align(fq,
          genomeDir,
          gtf,
          samplename,
          outdir,
          region,
          sc5p,
          core=4,
          logger=None,
          star_path='STAR',
          samtools_path='samtools',
          **kwargs):
    if ('steps' not in kwargs) or (not kwargs['steps']):
        kwargs['steps'] = ['STAR', 'SortByPos', 'FeatureCounts', 'SortByName']

    globals()['logger'] = logger

    basedir = os.path.join(outdir, 'step2')
    STAR_dir = os.path.join(basedir, 'STAR')
    os.makedirs(STAR_dir, exist_ok=True)
    prefix = os.path.join(STAR_dir, samplename + '_')

    if 'STAR' not in kwargs['steps']:
        logger.info('STAR skiped!')
    else:
        logger.info('STAR started!')
        bam, STARLog = run_STAR(fq=fq,
                                core=core,
                                genomeDir=genomeDir,
                                prefix=prefix,
                                star_path=star_path)
        logger.info('STAR done!')
    bam, STARLog =f'{prefix}Aligned.out.bam', f'{prefix}Log.final.out'
    # sort by pos
    if 'SortByPos' not in kwargs['steps']:
        logger.info('SortByPos skiped!')
    else:
        logger.info('SortByPos started!')
        bam = run_bamsort(
                    bam,
                    f'{prefix}SortedByCoordinate.bam',
                    core=core,
                    samtools_path=samtools_path)
        logger.info('SortByPos done!')
    bam = f'{prefix}SortedByCoordinate.bam'

    logger.info('run_qualimap started!')
    RnaSeqMetrics = run_qualimap(bam=bam, gtf=gtf, outdir=STAR_dir, SC5P=sc5p)
    logger.info('run_qualimap done!')

    with open(os.path.join(outdir, f'{samplename}_summary.json')) as fh:
        refpath=os.path.dirname(genomeDir.rstrip("/"))
        reffile=os.path.join(refpath,'reference.json')
        if os.path.exists(reffile):
                with open(reffile) as refjson:
                    refj=json.load(refjson)
                    genome=refj['genomes'][0]
        else:
                genome=genomeDir    
        summary = json.load(fh)
        summary['reference'] = genome
        Total = summary['stat']['total']
    summary_tmp = defaultdict()
    tmp = mapping_summary(STARLog, RnaSeqMetrics)
    Total = tmp['Number of input reads']
    mapped_genome_ratio = tmp['reads aligned']/Total
    summary_tmp['Reads Mapped to Genome'] = mapped_genome_ratio

    mapped_confident_ratio = (tmp['aligned to genes'] + tmp['no feature assigned']) / Total
    summary_tmp['Reads Mapped Confidently to Genome'] = mapped_confident_ratio

    mapped_intergenic_ratio = tmp['intergenic'] / Total
    summary_tmp['Reads Mapped to Intergenic Regions'] = mapped_intergenic_ratio
    
    mapped_intronic_ratio = tmp['intronic'] / Total
    summary_tmp['Reads Mapped to Intronic Regions'] = mapped_intronic_ratio

    mapped_exonic_ratio = tmp['exonic'] / Total
    summary_tmp['Reads Mapped to Exonic Regions'] = mapped_exonic_ratio

    with open(os.path.join(outdir, f'{samplename}_summary.json'), 'w') as fh:
        summary['mapping'] = summary_tmp
        json.dump(summary, fh, indent=4)


    featureCounts_dir = os.path.join(basedir, 'featureCounts')
    # FeatureCounts
    if 'FeatureCounts' not in kwargs['steps']:
        logger.info('run_featureCounts done!')
    else:
        os.makedirs(featureCounts_dir, exist_ok=True)
        logger.info('run_featureCounts started!')
        bam = run_featureCounts(bam=bam,
                                samplename=samplename,
                                outdir=featureCounts_dir,
                                gtf=gtf,
                                region=region,
                                SC5P=sc5p,
                                core=core)
        logger.info('run_featureCounts done!')
    bam = os.path.join(featureCounts_dir, f'{samplename}_SortedByCoordinate.bam.featureCounts.bam')
    # sort by name
    if 'SortByName' not in kwargs['steps']:
        logger.info('SortByName skiped!')
    else:
        logger.info('SortByName started!')
        bam = run_bamsort(
                    bam,
                    os.path.join(featureCounts_dir, f'{samplename}_SortedByName.bam'),
                    core=core,
                    byname=True,
                    samtools_path=samtools_path)
        logger.info('SortByName done!')
