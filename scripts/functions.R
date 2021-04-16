## FUNCTIONS

add_metadata <- function(data) {
  mito.genes <-
    grep(pattern = "^Mt\\.|^MT\\.|^mt\\.|^Mt-|^MT-|^mt-",
         x = rownames(x = GetAssayData(object = data)),
         value = TRUE)
  percent.mito <-
    Matrix::colSums(GetAssayData(object = data, slot = "counts")[mito.genes, ]) /
    Matrix::colSums(GetAssayData(object = data, slot = "counts"))
  data[['percent.mito']] <- percent.mito
  data[['percent.mito_log10']] <- log10(data[['percent.mito']] + 1)
  data[['nCount_RNA_log10']] <- log10(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log10']] <- log10(data[['nFeature_RNA']] + 1)
  data[['nCount_RNA_log2']] <- log2(data[['nCount_RNA']] + 1)
  data[['nFeature_RNA_log2']] <- log2(data[['nFeature_RNA']] + 1)
  data[['scaled_mito']] <- scale(percent.mito)
  data[['scaled_nCount_RNA']] <- scale(data[['nCount_RNA_log10']])
  attr(data$scaled_nCount_RNA, "scaled:center") <- NULL
  attr(data$scaled_nCount_RNA, "scaled:scale") <- NULL
  attr(data$scaled_mito, "scaled:center") <- NULL
  attr(data$scaled_mito, "scaled:scale") <- NULL
  data
}

draw_plots <- function(path, data) {
  VlnPlot(
    data,
    features = "nFeature_RNA",
    pt.size = 0.1
  ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_features.pdf'))
  VlnPlot(
    data,
    features = "nCount_RNA",
    pt.size = 0.1
  ) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'vln_plot_umi.pdf'))
  FeatureScatter(data, feature1 = "nCount_RNA_log10", feature2 = "nFeature_RNA_log10") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_features_log10_plot.pdf'))
  FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'umi_features_plot.pdf'))
  ElbowPlot(data, ndims = 50) + theme(aspect.ratio = 1)
  ggsave(paste0(path, 'elbow_plot.pdf'))
}

filter_mito <- function(obj, path, threshold=0.2) {
  VlnPlot(obj, features = "percent.mito")+
    theme(aspect.ratio = 1, legend.position = "none")+
    ggtitle(unique(obj$genotype))
  ggsave(paste0(paste0(path, unique(obj$genotype)), '_mt_before_filtration.pdf'))
  expr <- FetchData(object = obj, vars = 'percent.mito')
  obj <- obj[, which(x = expr < threshold)]
  VlnPlot(obj, features = "percent.mito")+
    theme(aspect.ratio = 1, legend.position = "none")+
    ggtitle(unique(obj$genotype))
  ggsave(paste0(paste0(path, unique(obj$genotype)), '_mt_after_filtration.pdf'))
  obj
}

analyze_object <- function(object, ident) {
  Idents(object) <- object[[ident]]
  cluster.averages <- AverageExpression(object = object, assays = 'SCT', slot = 'data')
  sapply(names(cluster.averages),
         function(x) write.table(cluster.averages[[x]], file=paste(x, paste(ident, "_clusters.tsv", sep = '_'), sep = '_')))

  whole.markers <- FindAllMarkers(object = object,
                                  assay='SCT',
                                  only.pos = TRUE,
                                  min.pct = 0.10,
                                  test.use = 'MAST')
  write.table(whole.markers, paste(ident, "markers.tsv", sep = '_'), sep="\t", quote=F, row.names=F)
}

get_df <- function(path) {
  sample_name <- gsub('_BAT.*', '', path)
  data <- Read10X(data.dir = path, gene.column = 2, unique.features = TRUE)
  bc_rank <- barcodeRanks(data)
  tot_counts <- colSums(data)
  tot_genes <- rowSums(data)
  data <- data[, tot_counts > metadata(bc_rank)$inflection]
  data <- data[tot_genes > 0, ]
  data <-
    CreateSeuratObject(
      counts = data,
      min.cells = 2,
      project = sample_name,
      min.features = 200
    )
  data <- add_metadata(data)
  data$sample <- sample_name
  data$genotype <- ifelse(grepl('WT', path), 'WT', 'KO')
  data
}

get_whole_obj <- function(pathes) {
  df <- lapply(pathes, function(x) get_df(x))
  names(df) <- sapply(df, function(x) unique(x$sample))
  c(df)
}
