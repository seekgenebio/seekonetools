suppressPackageStartupMessages({
    library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--indir", required = TRUE,
    help = "Required. The directory '*_feature_bc_matrix'.")
parser$add_argument("--name", required = TRUE,
    help = "Required. Sample name.")
parser$add_argument("--mincells", default = 0,
    help = "min.cells [default %(default)s]")
parser$add_argument("--minfeatures", default = 0,
    help = "min.features [default %(default)s]")
parser$add_argument("--outdir", default = "./",
    help = "The output directory [default %(default)s]")
parser$add_argument("--dims", default = 15,
    help = "PCA dims to use [default %(default)s]")
parser$add_argument("--minpct", default = 0.1,
    help = "min.pct  [default %(default)s]")
parser$add_argument("--logfc", default = 1,
    help = "logfc.threshold, base 2. [default %(default)s]")

Args <- parser$parse_args()

suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(ggplot2)
    library(jsonlite)
})

dir.create(Args$outdir)
setwd(Args$outdir)

data.count <- Read10X(data.dir = Args$indir)
features_df <-  read.table(file.path(Args$indir, 'features.tsv.gz'), sep="\t")
names(features_df)[1:2] <- c('Ensembl', 'gene')
data.dim <- dim(data.count)
obj <- CreateSeuratObject(counts = data.count, project = Args$name, min.cells = Args$mincells)
obj[['percent.mito']] <- PercentageFeatureSet(object = obj, pattern = '^MT-')
ggsave('VlnPlot.png', VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), pt.size=0.1 ,ncol = 3))

# summary_file <- normalizePath(sprintf('%s/../%s_summary.json', Args$outdir, Args$name))
# summary_json <- fromJSON(summary_file)
# summary_json$median_gene <- median(obj$nFeature_RNA)
# summary_json$median_umi <- median(obj$nCount_RNA)
# summary_json$median_mito <- median(obj$percent.mito)
# write(toJSON(summary_json, pretty=T, auto_unbox=T), summary_file)

plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave('FeatureScatter.png', CombinePlots(plots = list(plot1, plot2)))

write.table(round(do.call("cbind", tapply(obj$percent.mito, Idents(obj), quantile, probs=seq(0,1,0.05))), digits = 3),
            file='mito_quantile.xls', sep="\t", col.names=FALSE, quote=FALSE)
write.table(do.call("cbind", tapply(obj$nFeature_RNA, Idents(obj), quantile, probs=seq(0,1,0.05))),
            file='nFeature_quantile.xls', sep="\t", col.names=FALSE, quote=FALSE)
write.table(do.call("cbind", tapply(obj$nCount_RNA, Idents(obj), quantile,probs=seq(0,1,0.05))),
            file='nCount_quantile.xls', sep="\t", col.names=FALSE, quote=FALSE)

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

obj <- ScaleData(obj, vars.to.regress = "percent.mito")
obj <- RunPCA(obj)

# obj <- JackStraw(obj, num.replicate = 100)
# obj <- ScoreJackStraw(obj, dims = 1:20)

# ggsave('JackStrawPlot.png', JackStrawPlot(obj, dims = 1:15))
# ggsave('ElbowPlot.png', ElbowPlot(obj))

dims <- 1:Args$dims
obj <- FindNeighbors(obj, dims = dims)
obj <- FindClusters(obj, resolution = seq(0.2, 1.4, 0.3))

res.table <- sapply(grep('res', colnames(obj@meta.data), value = TRUE), function(x) length(unique(obj@meta.data[,x])))
write.table(res.table, file='resolution.xls', row.names=TRUE, col.names=FALSE, sep="\t", quote=FALSE)

Idents(obj) <- 'RNA_snn_res.0.8'
obj <- RunUMAP(obj, dims = dims)
ggsave('umap.png', DimPlot(obj, reduction = "umap"))

#obj <- RunTSNE(obj, dims = dims)
obj <- RunTSNE(obj, dims = dims,check_duplicates = FALSE)
ggsave('tsne.png', DimPlot(obj, reduction = "tsne"))

pal <-colorRampPalette(c("blue","cyan", "yellow","red"))
tsne_loci <- as.data.frame(Embeddings(obj, reduction='tsne'))
tsne_loci <- cbind(tsne_loci, obj[[]])

p <- ggplot(tsne_loci, aes(tSNE_1, tSNE_2, color=nCount_RNA))
p1 <- p + geom_point() + scale_colour_gradientn(colors=pal(500)) + theme_classic()
ggsave('tsne_umi.png', plot=p1)
write.table(tsne_loci, file='tsne_umi.xls', row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

obj.markers <- FindAllMarkers(obj, min.pct = Args$minpct, logfc.threshold = Args$logfc, only.pos = FALSE) %>%
    left_join(features_df, by='gene') %>% relocate(Ensembl, gene)
write.table(obj.markers, file='FindAllMarkers.xls', row.names=FALSE, sep="\t", quote=FALSE)

top10 <- obj.markers %>% group_by(cluster) %>% arrange(p_val, desc(avg_log2FC), desc(pct.1)) %>% head(n = 10)
obj <- ScaleData(obj, features=top10$gene)
png('top10_heatmap.png')
DoHeatmap(obj, features = top10$gene) + NoLegend()
dev.off()
