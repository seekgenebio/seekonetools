import os
from subprocess import run

def do_seurat(matrix,
              samplename,
              outdir,
              dims=15,
              minpct=0.1,
              logfc=1,
              logger=None,
              **kwargs):
    logger.info('seurat started!')
    outdir1 = os.path.join(outdir, 'step4')
    os.makedirs(outdir1, exist_ok=True)

    Rapp = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        'do_seurat.R')

    args = [
        'Rscript', Rapp, '--indir', matrix, '--name', samplename, '--outdir',
        outdir1, '--dims', dims, '--minpct', minpct, '--logfc', logfc
    ]
    args = [str(_) for _ in args]
    run(args, check=False)

    logger.info('seurat done!')
