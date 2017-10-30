library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))

threads <- 2

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- TRUE
} else {
    do_parallel <- FALSE
}

load("./cache/xpr.RData")
load("./cache/cn_genes.RData")
load("./cache/ne_genes.RData")
load("./cache/cce_genes.RData")

cache_dir <- "./data/ceres_cache_final"
dir.create(cache_dir, recursive=T, showWarnings=F)


ceres_input_wang <- load_ceres_inputs(
    inputs_dir="./data/ceres_inputs/wang/",
    replicates=F)
ceres_output_wang <- load_ceres_outputs(
    output_dir="./data/ceres_outputs/wang/")

wang_gene_level <- create_gene_level_df(ceres_input_wang, ceres_output_wang,
                                        cn_genes, xpr,
                                        ne_genes, cce_genes,
                                        do_parallel=do_parallel)

saveRDS(ceres_input_wang,
        file.path(cache_dir, "ceres_input_wang.rds"))
saveRDS(ceres_output_wang,
        file.path(cache_dir, "ceres_output_wang.rds"))
saveRDS(wang_gene_level,
        file.path(cache_dir, "wang_gene_level.rds"))

ceres_output_wang$guide_activity %>%
    write_tsv(file.path(cache_dir, "wang_guide_activity.tsv"))
wang_gene_level %>%
    select(Gene, CellLine, Gene_Effect) %>%
    df.to.mat() %>%
    write.csv(file.path(cache_dir, "wang_ceres_gene_effects.csv"))



ceres_input_gecko <- load_ceres_inputs(
    inputs_dir="./data/ceres_inputs/gecko/",
    replicates=F)
ceres_output_gecko <- load_ceres_outputs(
    output_dir="./data/ceres_outputs/gecko/")


gecko_gene_level <- create_gene_level_df(ceres_input_gecko, ceres_output_gecko,
                                        cn_genes, xpr,
                                        ne_genes, cce_genes,
                                        do_parallel=do_parallel)

saveRDS(ceres_input_gecko,
        file.path(cache_dir, "ceres_input_gecko.rds"))
saveRDS(ceres_output_gecko,
        file.path(cache_dir, "ceres_output_gecko.rds"))
saveRDS(gecko_gene_level,
        file.path(cache_dir, "gecko_gene_level.rds"))

ceres_output_gecko$guide_activity %>%
    write_tsv(file.path(cache_dir, "gecko_guide_activity.tsv"))
gecko_gene_level %>%
    select(Gene, CellLine, Gene_Effect) %>%
    df.to.mat() %>%
    write.csv(file.path(cache_dir, "gecko_ceres_gene_effects.csv"))




ceres_input_avana <- load_ceres_inputs(
    inputs_dir="./data/ceres_inputs/avana/",
    replicates=F)
ceres_output_avana <- load_ceres_outputs(
    output_dir="./data/ceres_outputs/avana/")


avana_gene_level <- create_gene_level_df(ceres_input_avana, ceres_output_avana,
                                        cn_genes, xpr,
                                        ne_genes, cce_genes,
                                        do_parallel=do_parallel)

saveRDS(ceres_input_avana,
        file.path(cache_dir, "ceres_input_avana.rds"))
saveRDS(ceres_output_avana,
        file.path(cache_dir, "ceres_output_avana.rds"))
saveRDS(avana_gene_level,
        file.path(cache_dir, "avana_gene_level.rds"))

ceres_output_avana$guide_activity %>%
    write_tsv(file.path(cache_dir, "avana_guide_activity.tsv"))
avana_gene_level %>%
    select(Gene, CellLine, Gene_Effect) %>%
    df.to.mat() %>%
    write.csv(file.path(cache_dir, "avana_ceres_gene_effects.csv"))





