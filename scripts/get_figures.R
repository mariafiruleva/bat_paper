suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(plyr))
suppressMessages(library(MASS))
suppressMessages(library(argparse))
suppressMessages(library(dyno))

source('scripts/funs_figures.R')

parser <-
  ArgumentParser(description = 'Get scRNA-seq related figures from the paper')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to seurat rda object')
parser$add_argument('--traj',
                    type = "character",
                    help = 'Path to trajectory rda object')
parser$add_argument('--pws',
                    type = "character",
                    help = 'Path to pathways file in json format')
parser$add_argument('--target_genes',
                    type = "character",
                    help = 'Path to txt file with gene names')
## SET VARIABLES

arguments <- parser$parse_args()

print(arguments)

load(arguments$data)
load(arguments$traj)

getPalette.1 <- colorRampPalette(brewer.pal(9, "Set1"))

## UMAP: clusters, split by genotype

whole.integrated$custom_clusters <- whole.integrated$integrated_snn_res.1

whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '0',
                                           '0', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '11',
                                           '1', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '2',
                                           '2', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '9',
                                           '3', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '1',
                                           '4', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '8',
                                           '5', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '7',
                                           '6', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '12',
                                           '7', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '13',
                                           '8', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '10',
                                           '9', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '5' |
                                             whole.integrated$integrated_snn_res.1 == '6' |
                                             whole.integrated$integrated_snn_res.1 == '14',
                                           '10', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '4' |
                                             whole.integrated$integrated_snn_res.1 == '16',
                                           '11', whole.integrated$custom_clusters)
whole.integrated$custom_clusters <- ifelse(whole.integrated$integrated_snn_res.1 == '3' |
                                             whole.integrated$integrated_snn_res.1 == '15',
                                           '12', whole.integrated$custom_clusters)

sorted_labels <- 0:12
whole.integrated$custom_clusters <- factor(x = whole.integrated$custom_clusters, levels = sorted_labels)
Idents(whole.integrated) <- 'custom_clusters'
new_labels <- c('0 Neutrophils', expression(~1~Myeloid~cells~italic(Zeb2)^{"hi"}),
                expression(~2~Monocytes~italic(Ly6C)^{"low"}), expression(~3~Monocytes~italic(Ly6C)^{"int"}),
                expression(~4~Monocytes~italic(Ly6C)^{"hi"}), '5 Matrix Macrophages',
                '6 Macrophages M2-like', expression(7~Macrophages~italic(Lpl)^{"hi"}),
                expression(8~Macrophages~italic(Plin2)^{"hi"}), '9 Dendritic cells',
                '10 T cells', '11 B cells', '12 NK cells')
whole.integrated$genotype <- factor(x = whole.integrated$genotype, levels = c("WT", "KO"))
plt <- DimPlot(whole.integrated, split.by='genotype', pt.size=0.25)+
  scale_color_manual(values=getPalette.1(length(unique(whole.integrated$custom_clusters))), labels = new_labels)+
  theme_bw(base_size=11)+
  theme(legend.text.align = 0, legend.key.size=unit(0.2, "in"), aspect.ratio = 1,
        plot.margin=grid::unit(c(0,0,0.2,0), "in"),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill = guide_legend(ncol = 4, override.aes = list(color = c('black'))))+
  scale_fill_continuous(guide="legend",breaks=seq(0.2,0.8,by=0.1))

LabelClusters(plt, id = "ident", size=5, fontface = 'bold', repel=F)


ggsave("article_plots/clustering_total.png", width = 8, height = 4, dpi=600, units='in')
ggsave("article_plots/clustering_total.svg", width = 11, height = 5, dpi=600, units='in', device = 'svg')

## histogram: cells per cluster per genotype

df <- cbind(as.character(whole.integrated$genotype), as.character(whole.integrated$custom_clusters)) %>%
  as.data.frame() %>%
  magrittr::set_colnames(c('genotype', 'cluster'))

df$cluster <- factor(x = df$cluster, levels = as.character(sort(as.numeric(levels(df$cluster)))))
df$genotype <- factor(x = df$genotype, levels = c("WT", "KO"))
p <- df %>% group_by(cluster, genotype) %>% dplyr::summarise(n = n())
per_gt <- table(whole.integrated$genotype) %>% as.data.frame() %>% magrittr::set_colnames(c('genotype', 'count'))

p <- p %>%
  group_by(genotype) %>%
  mutate(total = per_gt$count[match(genotype, per_gt$genotype)]) %>%
  group_by(cluster, add=TRUE) %>%
  mutate(per=round(100*n/total,2))

hist_plot <- ggplot(p,aes(x=cluster,y=per, fill=genotype))+
  geom_bar(stat = 'identity', position=position_dodge(width = 0.5), width = 0.45) +
  theme_bw(base_size=11)+
  scale_fill_brewer(palette="Set1", direction = -1)+
  theme(legend.text.align = 0, legend.key.size=unit(0.2, "in"), aspect.ratio = 1, legend.position = 'top', legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "in"))+
  geom_segment(aes(x = 2, y = 8, xend = 2, yend = 6),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 6, y = 13, xend = 6, yend = 11),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 7, y = 13, xend = 7, yend = 11),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 8, y = 10, xend = 8, yend = 8),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  geom_segment(aes(x = 9, y = 10, xend = 9, yend = 8),
               arrow = arrow(length = unit(0.3, "cm")), lineend='round')+
  ylab('cells per clutser, %')


ggsave(plot = hist_plot, filename = "article_plots/corrected_hist.png", width = 4, height = 4, dpi=600, units='in')
ggsave(plot = hist_plot, filename = "article_plots/corrected_hist.svg", width = 4, height = 4, dpi=600, units='in')

## Violin plots

get_expr <- function(object, gene_set, slot='data', assay='SCT') {
  av <- numeric(ncol(object))
  zz <- which(tolower(rownames(GetAssayData(object, slot = slot, assay = assay))) %in% tolower(gene_set))
  object@assays$SCT@data[zz, ]
}

plot_vln <- function(object, gene_set, reduction="umap", assay='SCT', slot='data') {
  red <- cbind(get_expr(object, gene_set)) %>% as.data.frame() %>%
    magrittr::set_colnames(c('expression'))
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c("WT", "KO"))
  red$cluster <- object$custom_clusters
  target_clusters <- c(2, 4, 3, 6, 5, 7, 8)
  red <- red %>% dplyr::filter(cluster %in% target_clusters)
  ggplot(data=red, aes(x=cluster, y=expression, fill=genotype)) + theme_bw(base_size = 8) +
    theme(panel.spacing = unit(0, "lines"),
          legend.position = 'none',
          legend.title = element_blank(),
          aspect.ratio = 0.2,
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "italic", size=6, margin=margin(0,0,0,0)),
          plot.margin=grid::unit(c(0,0,0,0), "in")) +
    geom_split_violin(scale="width", alpha=0.7) +
    scale_fill_brewer(palette='Set1', guide=guide_legend(ncol=2), direction=-1) +
    ylab(NULL) + xlab(NULL)+ggtitle(gene_set)
  if (!dir.exists('article_plots/violins')) {
    dir.create('article_plots/violins', recursive = T)
  }
  ggsave(sprintf("article_plots/violins/%s_vln.png", gene_set), width = 3, height = 0.8, dpi=600, units='in')
}

sapply(c('Plin2', 'Lpl', 'Cd36', 'Trem2', 'Fabp4', 'Fabp5'), function(x) plot_vln(whole.integrated, x))

## Trajectory plot

obj$cluster$grouping <- as.factor(obj$cluster$grouping)
obj$cluster$grouping <- ifelse(obj$cluster$grouping == 1, '4',
                               ifelse(obj$cluster$grouping == 9, '3',
                                      ifelse(obj$cluster$grouping == 8, '5',
                                             ifelse(obj$cluster$grouping == 7, '6',
                                                    ifelse(obj$cluster$grouping == 12, '7',
                                                    ifelse(obj$cluster$grouping == 13, '8', '2'))))))
obj$cluster$grouping <- as.factor(obj$cluster$grouping)
new_traj_labels <- c(expression(~2~Monocytes~italic(Ly6C)^{"low"}), expression(~3~Monocytes~italic(Ly6C)^{"int"}),
                     expression(~4~Monocytes~italic(Ly6C)^{"hi"}), '5 Matrix Macrophages', '6 Macrophages M2-like',
                     expression(7~Macrophages~italic(Lpl)^{"hi"}), expression(8~Macrophages~italic(Plin2)^{"hi"}))

traj <- plot_dimred(
  model,
  dimred = dimred,
  expression_source = obj$cluster$expression,
  color_density = "grouping",
  alpha_cells = 0.7,
  size_cells = 1,
  size_trajectory = 1.5,
  grouping = obj$cluster$grouping,
  label_milestones = F,
)+
  scale_fill_manual(values=getPalette.1(length(unique(whole.integrated$custom_clusters)))[3:9], labels = NULL)+
  scale_color_manual(values=getPalette.1(length(unique(whole.integrated$custom_clusters)))[3:9],labels = new_traj_labels)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=3)))+
  theme_void(base_size = 11)+
  theme(aspect.ratio = 1, legend.title = element_blank(),
        legend.text.align = 0, legend.key.size=unit(0.2, "in"), plot.margin=grid::unit(c(0,0,0.2,0), "in"),
        legend.text=element_text(size=11))

ggsave(plot = traj, filename = "article_plots/trajectory.png", width = 8, height = 6, dpi=600, units='in')
ggsave(plot = traj, filename = "article_plots/trajectory.svg", width = 8, height = 6, dpi=600, units='in', device = 'svg')


## PATHWAY

get_pw_expr <- function(object, gene_set, slot='data', assay='SCT') {
  av <- numeric(ncol(object))
  zz <- which(tolower(rownames(GetAssayData(object, slot = slot, assay = assay))) %in% tolower(gene_set))
  geneExp <- as.matrix(log2(object@assays$SCT@data[zz, ] + 1))
  geneExp <- t(scale(t(geneExp)))
  geneExp[is.nan(geneExp)] <- 0
  av <- av + colSums(geneExp) / length(gene_set)
  av
}


plot_target_pw <- function(object, gene_set, pw_name, reduction="umap", assay='SCT', slot='data', macs=F) {
  red <- cbind(get_pw_expr(object, gene_set), object@reductions[['umap']]@cell.embeddings) %>% as.data.frame() %>%
    magrittr::set_colnames(c('expression', paste0(reduction, 1), paste0(reduction, 2)))
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c('WT', 'KO'))
  p <- ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.15)+
    theme_bw(base_size=8) +
    facet_grid(.~genotype)+
    scale_color_gradientn(colours=c("darkblue", "blue", "grey", "red", "darkred"),
                          breaks = c(0, floor(max(red$expression) * 100) / 100),
                          rescaler = function(x, from) {
      res <- numeric(length(x))
      res[x >= 1 & !is.na(x)] <- 1
      res[x < 1 & !is.na(x)] <- (x[x < 1 & !is.na(x)] + 1) / 2
      res[x < -1 & !is.na(x)] <- 0
      res[is.na(x)] <- NA
      res
    })+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5),
          plot.margin=grid::unit(c(0,0,0,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gsub('_', ' ', pw_name))
  if (macs) {
    p <- p+
      xlim(NA, 5)+
      ylim(0, NA)
  }
  if (!dir.exists('article_plots/macs_pw')) {
    dir.create('article_plots/macs_pw', recursive = T)
  }
  ggsave(plot=p, filename=sprintf("article_plots/macs_pw/%s.png", pw_name), width = 4.25, height = 2, dpi=600, units='in')
}


pws <- jsonlite::fromJSON(arguments$pws)

## zoom macs: pathways

plot_target_pw(whole.integrated, pws$KEGG_GLYCEROLIPID_METABOLISM, 'KEGG_GLYCEROLIPID_METABOLISM', macs = T)
plot_target_pw(whole.integrated, pws$KEGG_GLYCEROPHOSPHOLIPID_METABOLISM, 'KEGG_GLYCEROPHOSPHOLIPID_METABOLISM', macs = T)
plot_target_pw(whole.integrated, pws$PID_LYSOPHOSPHOLIPID_PATHWAY, 'PID_LYSOPHOSPHOLIPID_PATHWAY', macs = T)
plot_target_pw(whole.integrated, pws$HALLMARK_FATTY_ACID_METABOLISM, 'HALLMARK_FATTY_ACID_METABOLISM', macs = T)

## zoom macs : genes

plot_target_gene <- function(object, gene, reduction="umap", assay='SCT', slot='data', macs = F) {
  data <- GetAssayData(object, slot = slot, assay = assay)
  red <- object@reductions[[reduction]]@cell.embeddings
  red <- as.data.frame(red)
  colnames(red) <- paste0(reduction, 1:ncol(red))
  genes_indexes <- which(rownames(data) %in% gene[[1]])
  expression_signaling <- data[genes_indexes,]
  red$expression <- expression_signaling
  red$genotype <- ifelse(rownames(red) %in% names(object$genotype[object$genotype == 'WT']), 'WT', 'KO')
  red$genotype <- factor(x = red$genotype, levels = c('WT', 'KO'))
  red$cluster <- object$integrated_snn_res.0.2
  onlyBreak <- floor(max(red$expression))
  p <- ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.1)+
    theme_bw(base_size=11) +
    facet_grid(.~genotype)+
    scale_color_gradientn(colours=c("grey", "red", "red3"),breaks=c(0, onlyBreak),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "italic"),
          plot.margin=grid::unit(c(0,0,0.2,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gene)
  if (macs) {
    p <- p+
      xlim(NA, 5)+
      ylim(0, NA)
  }
  if (!dir.exists('article_plots/macs_genes')) {
    dir.create('article_plots/macs_genes', recursive = T)
  }
  ggsave(sprintf("article_plots/macs_genes/%s.png", gene), width = 4.25, height = 2.2, dpi=600, units='in')
}

sapply(c('Mrc1', 'Clec10a', 'Mki67', 'Ccna2', 'Top2a',
         'Cd86', 'Cd80', 'Cd68', 'Tlr2', 'Tlr4', 'Cxcl10'),
       function(x) plot_target_gene(whole.integrated, x, macs = T))

## all genes

plot_markers <- function(object, gene, out_dir, reduction="umap", assay='SCT', slot='data') {
  data <- GetAssayData(object, slot = slot, assay = assay)
  red <- object@reductions[[reduction]]@cell.embeddings
  red <- as.data.frame(red)
  colnames(red) <- paste0(reduction, 1:ncol(red))
  genes_indexes <- which(rownames(data) %in% gene[[1]])
  expression_signaling <- data[genes_indexes,]
  red$expression <- expression_signaling
  onlyBreak <- floor(max(red$expression))
  ggplot(red, aes_string(x=paste0(reduction, 1),
                         y=paste0(reduction, 2),
                         color="expression"))+
    geom_point(size=0.05)+
    theme_bw(base_size=11) +
    scale_color_gradientn(colours=c("grey", "red", "red3"),breaks=c(0, onlyBreak),
                          guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
    theme(legend.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_blank(),
          aspect.ratio = 1,
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "italic"),
          plot.margin=grid::unit(c(0,0,0.2,0), "in"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank())+
    ggtitle(gene)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }
  ggsave(sprintf("%s/%s.png", out_dir, gene), width = 4.25, height = 2.2, dpi=600, units='in')
}

top5_markers <- data.table::fread(arguments$target_genes, header = F)

sapply(top5_markers$V1, function(x) plot_markers(whole.integrated, x, 'article_plots/marker_genes'))
