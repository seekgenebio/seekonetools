import os
import json
import click
from .chemistry import CHEMISTRY


_steps = {
    'step1': [],
    'step2': ['STAR', 'SortByPos', 'FeatureCounts', 'SortByName'],
    'step3': [],
    'step4': []
}

@click.group(help="quantifies 3' single-cell gene expression")
@click.option('--steps', default=None, type=click.Path(), help='')
@click.pass_obj
def rna(obj, steps):
    if  steps:
        from yaml import load
        try:
            from yaml import CLoader as Loader
        except ImportError:
            from yaml import Loader
        with open(steps) as fh:
            obj['steps'] = load(fh, Loader=Loader)
    else:
        obj['steps'] = _steps


@rna.command(help="extract cell barcode and umi.")
@click.option('--fq1', 'fq1', required=True, type=click.Path(), multiple=True, help='read1 fq file, can specify multiple times.')
@click.option('--fq2', 'fq2', required=True, type=click.Path(), multiple=True, help='read2 fq file, can specify multiple times.')
@click.option('--samplename', required=True, help='sample name.')
@click.option('--outdir', default='./', show_default=True, type=click.Path(), help='output dir.')
@click.option('--shift', is_flag=True, default=False, show_default=True, help='shift')
@click.option('--pattern', 'shift_pattern', default='A', help='')
@click.option('--barcode', multiple=True, help='barcode white list file, can specify multiple times.')
@click.option('--structure', help='')
@click.option('--linker', multiple=True, help='linker white list file, can specify multiple times.')
@click.option('--misB', 'B', nargs=2, default=(1, 0), type=click.Tuple([int, int]), show_default=True, help='err and indel')
@click.option('--misL', 'L', nargs=2, default=(1, 0), type=click.Tuple([int, int]), show_default=True, help='err and indel')
@click.option('--core', default=4, show_default=True, help='core')
@click.option('--chemistry', help='eg: SO01V3')
@click.pass_obj
def step1(obj, **kwargs):
    if kwargs['chemistry']:
        kwargs.update(CHEMISTRY[kwargs['chemistry']])

    from .step1 import barcode
    kwargs['logger'] = obj['logger']
    barcode(**kwargs)

@rna.command(help="align reads to genome.")
@click.option('--fq', required=True, help='')
@click.option('--genomeDir', 'genomeDir', required=True, type=click.Path(), help='')
@click.option('--gtf', required=True, type=click.Path(), help='')
@click.option('--samplename', required=True, help='')
@click.option('--outdir', default='./', show_default=True, type=click.Path(), help='')
@click.option('--star_path', 'star_path', default='STAR', help='')
@click.option('--samtools_path', 'samtools_path', default='samtools', help='')
@click.option('--core', default=4, show_default=True, help='')
@click.pass_obj
def step2(obj, **kwargs):
    from .step2 import align
    kwargs['logger'] = obj['logger']
    align(**kwargs)

# def count(bam, outdir, samplename, gtf, logger, expectNum=3000, **kwargs):
@rna.command(help="quantifies.")
@click.option('--bam', required=True, help='')
@click.option('--outdir', default='./', show_default=True, type=click.Path(), help='')
@click.option('--samplename', required=True, help='')
@click.option('--gtf', required=True, type=click.Path(), help='')
@click.option('--expectNum', 'expectNum', default=3000, show_default=True, help='')
@click.pass_obj
def step3(obj, **kwargs):
    from .step3 import count
    kwargs['logger'] = obj['logger']
    count(**kwargs)


@rna.command(help="seurat.")
@click.option('--matrix', type=click.Path(), help='')
@click.option('--samplename', required=True, help='')
@click.option('--outdir', default='./', show_default=True, type=click.Path(), help='')
@click.option('--dims', default=15, show_default=True, help='')
@click.option('--minpct', default=0.1, show_default=True, help='')
@click.option('--logfc', default=1.0, show_default=True, help='')
@click.option('--rscript_path', 'rscript_path', default='Rscript', help='')
@click.pass_obj
def step4(obj, **kwargs):
    from .step4 import do_seurat
    kwargs['logger'] = obj['logger']
    do_seurat(**kwargs)


@rna.command(help="report.")
@click.option('--samplename', required=True, help='')
@click.option('--outdir', default='./', show_default=True, type=click.Path(), help='')
@click.pass_obj
def report(obj, **kwargs):
    from .report import report
    kwargs['logger'] = obj['logger']
    report(**kwargs)


@rna.command(help="run all steps.")
@click.pass_obj
@click.option('--fq1', 'fq1', required=True, type=click.Path(), multiple=True, help='read1 fq file, can specify multiple times.')
@click.option('--fq2', 'fq2', required=True, type=click.Path(), multiple=True, help='read2 fq file, can specify multiple times.')
@click.option('--samplename', required=True, help='sample name.')
@click.option('--outdir', default='./', show_default=True, type=click.Path(), help='output dir.')
@click.option('--shift', is_flag=True, default=False, help='shift')
@click.option('--pattern', 'shift_pattern', default='A', help='')
@click.option('--barcode', multiple=True, help='barcode white list file, can specify multiple times.')
@click.option('--structure', help='')
@click.option('--linker', multiple=True, help='linker white list file, can specify multiple times.')
@click.option('--misB', 'B', nargs=2, default=(1, 0), type=click.Tuple([int, int]), show_default=True, help='')
@click.option('--misL', 'L', nargs=2, default=(1, 0), type=click.Tuple([int, int]), show_default=True, help='')
@click.option('--core', default=4, show_default=True, help='')
@click.option('--genomeDir', 'genomeDir', required=True, type=click.Path(), help='')
@click.option('--gtf', required=True, type=click.Path(), help='')
@click.option('--star_path', 'star_path', default='STAR', help='')
@click.option('--samtools_path', 'samtools_path', default='samtools', help='')
@click.option('--rscript_path', 'rscript_path', default='Rscript', help='')
@click.option('--chemistry', help='eg: SO01V3')
@click.option('--expectNum', 'expectNum', default=3000, show_default=True, help='')
def run(obj, **kwargs):
    if kwargs['chemistry']:
        kwargs.update(CHEMISTRY[kwargs['chemistry']])

    kwargs['logger'] = obj['logger']
    kwargs['outdir'] = os.path.join(kwargs['outdir'], kwargs['samplename'])
    if 'step1' in obj['steps']:
        from .step1 import barcode 
        barcode(**kwargs)
    fq = os.path.join(kwargs['outdir'], 'step1', f'{kwargs["samplename"]}_2.fq.gz')
    kwargs['fq'] = fq

    if 'step2' in obj['steps']:
        from .step2 import align
        align(flags=_steps['step2'], **kwargs)
    bam = os.path.join(kwargs['outdir'], 'step2', 'featureCounts',  f'{kwargs["samplename"]}_SortedByName.bam')
    kwargs['bam'] = bam

    if 'step3' in obj['steps']:
        from .step3 import count
        count(**kwargs)

    matrix = os.path.join(kwargs['outdir'], 'step3', 'filtered_feature_bc_matrix')
    kwargs['matrix'] = matrix
    if 'step4' in obj['steps']:
        from .step4 import do_seurat
        do_seurat(**kwargs)

        from .report import report
        report(**kwargs)

