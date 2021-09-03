import click

@click.group(help="utils.")
@click.pass_context
def utils(ctx):
    pass

@utils.command(help="gene type summary.")
@click.option('--gtf', type=click.Path(), help="gtf file.")
@click.option('--feature', default='gene', show_default=True, help="feature, e.g., gene, transcript")
@click.option('--key', default='gene_biotype', show_default=True, help='attribution key, e.g., gene_biotype')
def gtfstat(gtf, feature, key='gene_biotype'):
    from .mkref import gtfstat
    gtfstat(gtf, feature, key)  
    
@utils.command(help="gtf filter.")
@click.option('--gtf', type=click.Path(), help="gtf file.")
@click.option('--biotype',  multiple=True, help="biotype to keep.")
@click.option('--key', default='gene_biotype', show_default=True, help='attribution key, e.g., gene_biotype')
def gtffilter(gtf, biotype, key='gene_biotype'):
    from .mkref import gtffilter
    gtffilter(gtf, biotype, key)

@utils.command(help="addtag.")
def addtag():
    pass
