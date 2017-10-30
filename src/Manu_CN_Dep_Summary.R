library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/cn_dep_summary", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)



avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")
gecko_gene_level <- readRDS("./data/ceres_cache_final/gecko_gene_level.rds")
wang_gene_level <- readRDS("./data/ceres_cache_final/wang_gene_level.rds")


rpkm_thresh <- -1


cn_dep_boxplot(avana_gene_level) +
    ggsave(file.path(out_dir, "cn_dep_boxplot_avana.pdf"), width=8, height=4)


cn_dep_correlation(avana_gene_level, title="All Genes") +
    ggsave(file.path(out_dir, "cn_dep_cor_avana.pdf"),
           width=4, height=3)

avana_unxpr <- avana_gene_level %>%
    filter(!is.na(RPKM), RPKM < rpkm_thresh)

cn_dep_boxplot(avana_unxpr) +
    ggsave(file.path(out_dir, "cn_dep_boxplot_avana_unxpr.pdf"), width=8, height=4)

avana_unxpr %>%
    group_by(Gene) %>%
    filter(n() >= 3) %>%
    ungroup %>%
    cn_dep_correlation(title="Unexpressed Genes") +
    ggsave(file.path(out_dir, "cn_dep_cor_avana_unxpr.pdf"),
           width=4, height=3)



cn_dep_boxplot(gecko_gene_level) +
    ggsave(file.path(out_dir, "cn_dep_boxplot_gecko.pdf"), width=8, height=4)

cn_dep_correlation(gecko_gene_level, title="All Genes") +
    ggsave(file.path(out_dir, "cn_dep_cor_gecko.pdf"),
           width=4, height=3)

gecko_unxpr <- gecko_gene_level %>%
    filter(!is.na(RPKM), RPKM < rpkm_thresh)

cn_dep_boxplot(gecko_unxpr) +
    ggsave(file.path(out_dir, "cn_dep_boxplot_gecko_unxpr.pdf"), width=8, height=4)

gecko_unxpr %>%
    group_by(Gene) %>%
    filter(n() >= 3) %>%
    ungroup %>%
    cn_dep_correlation(title="Unexpressed Genes") +
    ggsave(file.path(out_dir, "cn_dep_cor_gecko_unxpr.pdf"),
           width=4, height=3)



cn_dep_boxplot(wang_gene_level) +
    ggsave(file.path(out_dir, "cn_dep_boxplot_wang.pdf"), width=8, height=4)

cn_dep_correlation(wang_gene_level, title="All Genes") +
    ggsave(file.path(out_dir, "cn_dep_cor_wang.pdf"),
           width=4, height=3)

wang_unxpr <- wang_gene_level %>%
    filter(!is.na(RPKM), RPKM < rpkm_thresh)

cn_dep_boxplot(wang_unxpr) +
    ggsave(file.path(out_dir, "cn_dep_boxplot_wang_unxpr.pdf"), width=8, height=4)

wang_unxpr %>%
    group_by(Gene) %>%
    filter(n() >= 3) %>%
    ungroup %>%
    cn_dep_correlation(title="Unexpressed Genes") +
    ggsave(file.path(out_dir, "cn_dep_cor_wang_unxpr.pdf"),
           width=4, height=3)

