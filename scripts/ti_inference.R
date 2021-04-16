suppressMessages(library(argparse))
suppressMessages(library(babelwhale))
suppressMessages(library(dyno))
suppressMessages(library(dplyr))
suppressMessages(library(ggrepel))
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))


parser <-
  ArgumentParser(description = 'Trajectory analysis using dyno wrapper')
parser$add_argument('--assay_data',
                    type = "character",
                    help = 'Assay name: RNA / SCT / integrated')
parser$add_argument('--assay_counts',
                    type = "character",
                    help = 'Assay name: RNA / SCT / integrated')
parser$add_argument('--data',
                    type = "character",
                    help = 'Path to the Seurat object in RData format')
parser$add_argument('--clusters',
                    type = "integer", nargs='+',
                    help = 'Clusters ids')
parser$add_argument('--ident',
                    type = "character",
                    help = 'Ident name, e.g., integrated_snn_res.0.6')
parser$add_argument('--out_dir',
                    type = "character",
                    help = 'path to output dir')
parser$add_argument('--ti_tool',
                    type = "character",
                    help = 'Name of the trajectory inference tool available at dyno wrapper')
parser$add_argument('--sng_cache',
                    type = "character",
                    help = 'Path to the cahce dir for singularity container')

## SET VARIABLES

args <- parser$parse_args()
clusters <- args$clusters
ti_tool <- args$ti_tool
setwd(args$out_dir)
print(args)

## CONFIG SINGULARITY

babelwhale::set_default_config(babelwhale::create_singularity_config(cache_dir = args$sng_cache))

## DEFINE SOME FUCNTIONS

get_data <- function(path, clusters, ident, assay_counts='SCT', assay_data='integrated') {
  obj <- get(load(path))
  Idents(obj) <- ident
  obj <- subset(obj, idents = clusters)
  wrap_data <- wrap_expression(
    counts = t(as.matrix(obj@assays[[assay_counts]]@counts)),
    expression = t(as.matrix(obj@assays[[assay_data]]@data))
  )
  wrap_data_by_clust <- add_grouping(
    wrap_data,
    obj[[ident]][[1]]
  )
  wrap_data_by_gt <- add_grouping(
    wrap_data,
    obj$genotype
  )
  res <- list(wrap_data_by_clust, wrap_data_by_gt)
  names(res) <- c('cluster', 'genotype')
  res
}

get_branching_point <- function(model, obj) {
  branching_milestone <- model$milestone_network %>% group_by(from) %>% filter(n() > 1) %>% pull(from) %>% dplyr::first()
  branch_point <- calculate_branching_point_feature_importance(model,
                                                               expression_source=obj$expression,
                                                               milestones_oi = branching_milestone)
  branch_point
}

get_hmap <- function(x, branch_feature_importance, model, wrap_data) {
  plot_heatmap(model, expression_source = obj$cluster$expression, features_oi =
                 branch_feature_importance %>%
                 filter(from == x[1]) %>%
                 filter(to == x[2]) %>%
                 top_n(50, importance) %>%
                 pull(feature_id) %>%
                 unique() %>% as.character())+
    ggtitle(sprintf('top-50, from %s to %s', x[1], x[2]))
  ggsave(sprintf('plots/hmaps/top50_from_%s_to-%s.png', x[1], x[2]))
  plot_heatmap(model, expression_source = obj$cluster$expression, features_oi =
                 branch_feature_importance %>%
                 filter(from == x[1]) %>%
                 filter(to == x[2]) %>%
                 top_n(100, importance) %>%
                 pull(feature_id) %>%
                 unique() %>% as.character())+
    theme(text = element_text(size=12))+
    ggtitle(sprintf('top-100, from %s to %s', x[1], x[2]))
  ggsave(sprintf('plots/hmaps/top100_from_%s_to-%s.png', x[1], x[2]),
         width=18, height=18)
}

save_tables <- function(model, branch_feature_importance, branch_point) {
  write.table(model$progressions,
              file='tables/progression.csv', quote = F, row.names = F, sep=',')
  write.table(model$milestone_network,
              file='tables/milestone_network.csv', quote = F, row.names = F, sep=',')
  write.table(branch_feature_importance,
              file='tables/branch_feature_importance.csv', quote = F, row.names = F, sep=',')
  write.table(branch_point,
              file='tables/branch_point.csv', quote = F, row.names = F, sep=',')
}

## PREPARE THE OBJECT

obj <- get_data(args$data, clusters, args$ident)

## INFER TRAJECTORY

set.seed(1)

model <- infer_trajectory(obj$cluster, method = args$ti_tool, verbose = TRUE)

## UMAP

dimred <- dyndimred::dimred_umap(obj$cluster$expression)

## CALCULATE SOME FEATURES: FEATURE IMPORTANCE, BRANCHING POINT

branch_feature_importance <- calculate_branch_feature_importance(model, expression_source=obj$cluster$expression)
branch_point <- get_branching_point(model, obj$cluster)

## VISUALIZATION

plot_dimred(
  model,
  dimred = dimred,
  expression_source = obj$cluster$expression,
  color_density = "grouping",
  alpha_cells = 0.7,
  grouping = obj$cluster$grouping,
  label_milestones = T,
)+
  theme(aspect.ratio = 1, legend.position = "none")+
  ggtitle(sprintf('%s, clusters: %s', args$ident, paste0(clusters, collapse=', ')))
ggsave('plots/trajectory_clusters.png')

plot_dimred(
  model,
  dimred = dimred,
  expression_source = obj$genotype$expression,
  color_density = "grouping",
  alpha_cells = 0.7,
  grouping = obj$genotype$grouping,
  label_milestones = T,
)+
  theme(aspect.ratio = 1, legend.position = "none")+
  ggtitle(sprintf('%s, %s, genotype: WT, KO', args$ident, paste0(clusters, collapse=', ')))
ggsave('plots/trajectory_genotype.png')

dir.create('plots/hmaps')
apply(model$milestone_network %>% select(from, to), 1, function(x) get_hmap(x, branch_feature_importance, model, obj$cluster))


## SAVE OUTPUT: TABLES

save_tables(model, branch_feature_importance, branch_point)

## SAVE OUTPUT: RDATA

save(list = c('args', 'obj', 'model', 'dimred', 'branch_feature_importance', 'branch_point'),
     file = 'rdata/trajectory.RData')

