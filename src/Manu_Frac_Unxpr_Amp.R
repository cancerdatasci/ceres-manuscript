library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/fraction_unxpr_amp", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")
gecko_gene_level <- readRDS("./data/ceres_cache_final/gecko_gene_level.rds")
wang_gene_level <- readRDS("./data/ceres_cache_final/wang_gene_level.rds")


avana_fps <- fraction_amps_unxpr(avana_gene_level, xlims=c(0,-8))

ggsave(file.path(out_dir, "fp_amplified_reldep_avana.pdf"),
       avana_fps$frac_amps_gg, width=4, height=4)

ggsave(file.path(out_dir, "fp_unexpressed_reldep_avana.pdf"),
       avana_fps$frac_unxpr_gg, width=4, height=4)


gecko_fps <- fraction_amps_unxpr(gecko_gene_level, xlims=c(0,-6))

ggsave(file.path(out_dir, "fp_amplified_reldep_gecko.pdf"),
       gecko_fps$frac_amps_gg, width=4, height=4)

ggsave(file.path(out_dir, "fp_unexpressed_reldep_gecko.pdf"),
       gecko_fps$frac_unxpr_gg, width=4, height=4)


wang_fps <- fraction_amps_unxpr(wang_gene_level, xlims=c(0,-6))

ggsave(file.path(out_dir, "fp_amplified_reldep_wang.pdf"),
       wang_fps$frac_amps_gg, width=4, height=4)

ggsave(file.path(out_dir, "fp_unexpressed_reldep_wang.pdf"),
       wang_fps$frac_unxpr_gg, width=4, height=4)


