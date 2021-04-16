import re

configfile: 'config.yaml'

rule all:
    input: 'article_plots/marker_genes'

rule seurat:
    output: rda=f'{config["init_path"]}/bat.RData'
    params: fast_tsne=config['fast_tsne'], rscript=config['r_script']
    log: 'logs/seurat.log'
    shell: '{params.rscript} scripts/analysis.R --fast_tsne {params.fast_tsne} 2> {log}'

rule ti_analysis:
    input: rda=rules.seurat.output.rda
    output: trajectory_rda="ti_inference/rdata/trajectory.RData",
          plt="ti_inference/plots/trajectory_clusters.png",
          tbls="ti_inference/tables/progression.csv",
          hmps=directory("ti_inference/plots/hmaps")
    params: sng_cache=config['sng_cache'], clst=lambda wildcards: re.sub(',', ' ', config['clusters_macs']),
            assay_counts=config['assay_counts'], assay_data=config['assay_data'], ident='integrated_snn_res.1',
            out_dir='ti_inference', ti_tool=config['ti_tool'], rscript=config['r_script']
    log: 'logs/ti_inference.log'
    shell:
         """
         {params.rscript} scripts/ti_inference.R --data {input.rda} \
         --assay_counts {params.assay_counts} --assay_data {params.assay_data} \
         --ident {params.ident} --ti_tool {params.ti_tool} --clusters {params.clst} \
         --sng_cache {params.sng_cache} --out_dir {params.out_dir} 2> {log}
         """

rule get_figures:
    input: rda=rules.seurat.output.rda, traj=rules.ti_analysis.output.trajectory_rda
    output: directory('article_plots/marker_genes')
    params: pws=config['pw_json'], rscript=config['r_script'], markers=config['markers_to_draw']
    log: 'logs/get_figures.log'
    shell: '{params.rscript} scripts/get_figures.R --data {input.rda} --traj {input.traj} --pws {params.pws} --genes {params.markers} 2> {log}'
