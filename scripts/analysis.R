suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(DropletUtils))
suppressMessages(library(functools))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(homologene))
suppressMessages(library(Matrix))
suppressMessages(library(sctransform))
suppressMessages(library(Seurat))
suppressMessages(library(argparse))

set.seed(1)

source('functions.R')

parser <-
  ArgumentParser(description = 'Seurat pipeline for BAT scRNA-seq project. Run in directory with filtered matrixes')
parser$add_argument('--fast_tsne',
                    type = "character",
                    help = 'Path to binary FIt-SNE package')
## SET VARIABLES

args <- parser$parse_args()

## CREATE DIRECTORY FOR PLOTS
path <- 'plots_stoyan/'
dir.create(path)

## GATHERING DATA TOGETHER

input_dirs <- list.dirs(full.names = F)[grepl("filtered_feature_bc_matrix", list.dirs())]

whole <- get_whole_obj(input_dirs)

## Number of cells before

cells.before <- sapply(whole, function(x) dim(GetAssayData(object = x, slot = "counts"))[2])

## FILTER MT CONTENT

whole <- sapply(whole, function(x) filter_mito(x, path))

## NORMALIZATION
whole <- sapply(whole, function(x) SCTransform(
  x,
  ncells=min(100000, ncol(x)),
  vars.to.regress = c("percent.mito"),
  verbose = T,
  conserve.memory = T
))

## INTEGRATION

whole.features <- SelectIntegrationFeatures(object.list = whole, nfeatures = 2000)
whole <- PrepSCTIntegration(object.list = whole, anchor.features = whole.features,
                            verbose = FALSE)
whole.anchors <- FindIntegrationAnchors(object.list = whole, normalization.method = "SCT",
                                        anchor.features = whole.features, verbose = FALSE)
whole.integrated <- IntegrateData(anchorset = whole.anchors, normalization.method = "SCT",
                                  verbose = FALSE)
gc()

## PCA

whole.integrated <- RunPCA(whole.integrated, verbose = FALSE)

## TSNE

whole.integrated <- RunTSNE(whole.integrated, dims = 1:20, tsne.method = "FIt-SNE",
                            fast_tsne_path = args$fast_tsne, nthreads = 4, max_iter = 2000)

## UMAP

whole.integrated <- RunUMAP(whole.integrated, dims = 1:20)

## CLUSTERING

whole.integrated <- FindNeighbors(object = whole.integrated, dims = 1:20)
whole.integrated <- FindClusters(object = whole.integrated, resolution = seq(0.2, 1, 0.2))

## Number of cells after

cells.after <- sapply(whole, function(x) length(colnames(x = x)))
cells.diff <- cells.before-cells.after
rbind(cells.before, cells.after, cells.diff)

## SAVING: RDATA

save(list = c('whole.integrated', 'whole.features', 'whole.anchors'), file = 'bat.RData')

## MARKERS

idents <- c('integrated_snn_res.0.2', 'integrated_snn_res.0.4',
            'integrated_snn_res.0.6', 'integrated_snn_res.0.8', 'integrated_snn_res.1')

sapply(idents, function(ident) analyze_object(object = whole.integrated, ident = ident))