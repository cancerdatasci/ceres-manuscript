library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/auc_sunburst", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)



load("./cache/cce_genes.RData")
load("./cache/ne_genes.RData")
load("./cache/lineage_df.RData")
avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")
gecko_gene_level <- readRDS("./data/ceres_cache_final/gecko_gene_level.rds")
wang_gene_level <- readRDS("./data/ceres_cache_final/wang_gene_level.rds")

avana_auc <- run_auc(avana_gene_level, lineage_df, ne_genes, cce_genes)
gecko_auc <- run_auc(gecko_gene_level, lineage_df, ne_genes, cce_genes)
wang_auc <- run_auc(wang_gene_level, lineage_df, ne_genes, cce_genes)


auc_sunburst(avana_auc, max_label_size=6, min_label_size=3) +
    ggsave(file.path(out_dir, "auc_sunburst_avana.pdf"), width=8, height=8)

auc_sunburst(gecko_auc, min_lines=1, min_label_size=5, max_label_size=6.5,
                  label_pos_y=0.85) +
    ggsave(file.path(out_dir, "auc_sunburst_gecko.pdf"), width=8, height=8)

auc_sunburst(wang_auc, min_lines=1, min_label_size=5, max_label_size=8,
                  label_pos_y=0.89, label_pos_x=0.2333333, fill_begin=0.4) +
    ggsave(file.path(out_dir, "auc_sunburst_wang.pdf"), width=8, height=8)
