import os
import base64
import click
import pandas as pd
import lzstring


template_dir = os.path.join(os.path.dirname(__file__), 'template')
data = {}


def diff_table(f, N=20):
    df = (pd.read_table(f)
         .assign(cluster=lambda df: df['cluster'].map(lambda x: f'cluster{x}'))
    )
    # table_id, classes
    fold_change_title = 'avg_logFC'
    if 'avg_log2FC' in df.columns:
        fold_change_title = 'avg_log2FC'
    df2 = (df.drop(columns=['p_val',  'pct.1',  'pct.2'])
        .groupby('cluster')
        .apply(lambda x: x.sort_values(['p_val_adj', fold_change_title], ascending=[True, False]).head(N))
        .reset_index(drop=True)
        .pivot_table(values=[fold_change_title, 'p_val_adj'], index=['Ensembl', 'gene'], columns='cluster')
        .swaplevel(axis=1)
        .sort_index(1, 0, key=lambda x: x.str.replace('cluster', '').astype(int))
    )
    # index_names
    return (df2.to_html(table_id='marker_table', classes='display', na_rep='-', index_names=False)
        .replace('border="1"', ''))


def barcode_rank_data(countsFile, barcodesFile):
    df = pd.read_csv(countsFile, sep='\t')
    barcodes = pd.read_csv(barcodesFile, header=None, sep='\t')
    UMIcounts = df.groupby(['cellID'])['UMINum'].agg(UMIcounts=sum)
    UMIcounts = UMIcounts.sort_values(by='UMIcounts', ascending=False)
    UMIcounts = UMIcounts.reset_index()
    #UMIcounts['rank'] = UMIcounts.index + 1

    #max_cellRank = int(UMIcounts.loc[UMIcounts['cellID'].isin(barcodes[0]), 'rank'].iloc[-1])
    max_idx_cell = int(UMIcounts[UMIcounts['cellID'].isin(barcodes[0])].index[-1])
    data['cells_index'] = [0, max_idx_cell]  # rows=5, index=[0,1,2,3,4], 3 cell, [0,2] => [0, 1, 2] is cell
    data['cells_data'] = UMIcounts.UMIcounts[0: max_idx_cell + 1].to_list() # [0: 2+1] => [0, 1, 2]

    # start from the last cell to prevent line gap
    data['background_index'] = [max_idx_cell, int(UMIcounts.index[-1])]  # [2, 4] => [2, 3, 4] is background, start from 2 (cell) to prevent line gap
    data['background_data'] = UMIcounts.UMIcounts[max_idx_cell:].to_list() # [2: ] => [2, 3, 4]

def reduction_data(reduction_umi):
    df = pd.read_table(reduction_umi, index_col=0)
    df = df.drop(['orig.ident', 'nFeature_RNA', 'percent.mito'], axis=1)
    data['reduction'] = df.to_dict(orient='list')
    data['params'] = [ _ for _ in df.columns if _.startswith(('RNA_snn_res','seurat_clusters')) ]
    data['labs'] = df.columns[0:2].to_list()

def png2base64(f):
    with open(f, 'rb') as fh:
        return base64.b64encode(fh.read()).decode()

def report(samplename, outdir, **kwargs):
    import json
    from jinja2 import Environment, FileSystemLoader

    summary_file = os.path.join(outdir, f'{samplename}_summary.json')
    assert os.path.exists(summary_file), f'{summary_file} not found!'
    with open(summary_file) as fh:
        summary = json.load(fh)

    sequencing_table = {}
    sequencing_table['Number of Reads'] = f'{summary["stat"]["total"]:,}'
    sequencing_table['Valid Barcodes'] = f'{summary["stat"]["valid"]/summary["stat"]["total"]:.2%}'
    sequencing_table['Sequencing Saturation'] = f'{summary["cells"]["Sequencing Saturation"]:.2%}'
    del summary["cells"]["Sequencing Saturation"]

    if 'no_anchor' in summary["stat"]:
        sequencing_table['Without Anchor'] = f'{summary["stat"]["no_anchor"]:,}'
    if 'B' in summary["stat"]:
        sequencing_table['Without Barcode'] = f'{summary["stat"]["B"]:,}'
    if 'L' in summary["stat"]:
        sequencing_table['Without Linker'] = f'{summary["stat"]["L"]:,}'
    if 'trimmed' in summary["stat"]:
        sequencing_table['Trimmed'] =  f'{summary["stat"]["trimmed"]:,}'
    if 'too_short' in summary["stat"]:
        sequencing_table['Too Short'] =  f'{summary["stat"]["too_short"]:,}'
    b_total_base = sum([sum(v) for v in summary["barcode_q"].values()])
    b30_base = sum([sum(v[30:]) for v in summary["barcode_q"].values()])
    sequencing_table['Q30 Bases in Barcode'] = f'{b30_base/b_total_base:.2%}'
    u_total_base = sum([sum(v) for v in summary["umi_q"].values()])
    u30_base = sum([sum(v[30:]) for v in summary["umi_q"].values()])
    sequencing_table['Q30 Bases in UMI'] = f'{u30_base/u_total_base:.2%}'

    mapping_table = {k: f'{v:.2%}' for k, v in summary["mapping"].items()}

    cells_table = dict([(k, f'{v:,}') if isinstance(v, int) else (k,f'{v:.2%}') for k,v in summary["cells"].items()])

    sample_table = {
        'Name': samplename, 
        'Description': '',
        'Transcriptome': summary["reference"],
        'Chemistry': '',
        'Seekone tools Version': summary["__version__"]
    }
    
    reduction_xls = os.path.join(outdir, 'step4', 'tsne_umi.xls')
    assert os.path.exists(reduction_xls), f'{reduction_xls} not found!'
    reduction_data(reduction_xls)

    count_xls = os.path.join(outdir, 'step3', 'counts.xls')
    barcodes_tsv = os.path.join(outdir, 'step3', 'filtered_feature_bc_matrix', 'barcodes.tsv.gz')
    barcode_rank_data(count_xls, barcodes_tsv)

    f =  os.path.join(outdir, 'step4', 'FindAllMarkers.xls')
    marker_table = diff_table(f)

    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template('base.html')
    with open(os.path.join(outdir, f'{samplename}_report.html'), 'w') as fh:
        rawdata = lzstring.LZString().compressToBase64(json.dumps(data))
        fh.write(template.render(
            sequencing_table = sequencing_table,
            mapping_table = mapping_table,
            cells_table = cells_table,
            sample_table = sample_table,
            marker_table = marker_table,
            downsample = summary["downsample"],
            rawdata=rawdata)
        )

    header=('Samplename,Estimated_Number_of_Cells,Mean_Reads_per_Cell,Median_Genes_per_Cell,Number_of_Reads,'
            'Valid_Barcodes,Sequencing_Saturation,Reads_Mapped_Confidently_to_Genome,Fraction_Reads_in_Cells,'
            'Total_Genes_Detected,Median_UMI_Counts_per_Cell')

    summary_data = [
             samplename,
             cells_table['Estimated Number of Cells'],
             sequencing_table['Number of Reads'],
             cells_table['Mean Reads per Cell'],
             cells_table['Median Genes per Cell'],
             sequencing_table['Number of Reads'],
             sequencing_table['Valid Barcodes'],
             sequencing_table['Sequencing Saturation'],
             mapping_table['Reads Mapped Confidently to Genome'],
             cells_table['Fraction Reads in Cells'],
             cells_table['Total Genes Detected'],
             cells_table['Median UMI Counts per Cell']
           ]

    with open(os.path.join(outdir, f'{samplename}_summary.csv'), 'w') as fh:
        fh.write(header + '\n')
        fh.write(','.join(str(_).replace(',', '') for _ in summary_data)+ '\n')

