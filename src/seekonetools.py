import click
from loguru import logger
from seekonetools.rna import rna
from seekonetools.utils import utils
from seekonetools.utils._version import __version__

@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version=__version__, message='%(version)s')
@click.pass_context
def cli(ctx):
    ctx.ensure_object(dict)
    ctx.obj['logger'] = logger
    ctx.obj['version'] = __version__

cli.add_command(rna)
cli.add_command(utils)

def main():
    cli()

# if __name__ == '__main__':
#     cli()
