library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))

threads <- 1

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- TRUE
} else {
    do_parallel <- FALSE
}

out_dir <- file.path("./output/figures/uncorrected_slopes", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


load("./cache/avana_neg_ctrl_guides.RData")
load("./cache/gecko_neg_ctrl_guides.RData")
load("./cache/wang_neg_ctrl_guides.RData")
load("./cache/avana_pos_ctrl_guides.RData")
load("./cache/gecko_pos_ctrl_guides.RData")
load("./cache/wang_pos_ctrl_guides.RData")
load("./cache/lineage_df.RData")
load("./cache/ccle_mut.RData")

ceres_input_avana <- readRDS("./data/ceres_cache_final/ceres_input_avana.rds")
ceres_input_gecko <- readRDS("./data/ceres_cache_final/ceres_input_gecko.rds")
ceres_input_wang <- readRDS("./data/ceres_cache_final/ceres_input_wang.rds")



avana_uncorrected_slopes <- run_uncorrected_slopes(ceres_input_avana, intercept=T,
                                             neg_ctrl_guides=avana_neg_ctrl_guides,
                                             pos_ctrl_guides=avana_pos_ctrl_guides,
                                             do_parallel=do_parallel,
                                             draw_plots=T,
                                             plot_dir=file.path(out_dir, "uncorrected_plots_avana"),
                                             dep_ylims=c(-8, 3))


gecko_uncorrected_slopes <- run_uncorrected_slopes(ceres_input_gecko, intercept=T,
                                                   neg_ctrl_guides=gecko_neg_ctrl_guides,
                                                   pos_ctrl_guides=gecko_pos_ctrl_guides,
                                                   do_parallel=do_parallel,
                                                   draw_plots=T,
                                                   plot_dir=file.path(out_dir, "uncorrected_plots_gecko"),
                                                   dep_ylims=c(-8, 3))

wang_uncorrected_slopes <- run_uncorrected_slopes(ceres_input_wang, intercept=T,
                                                  neg_ctrl_guides=wang_neg_ctrl_guides,
                                                  pos_ctrl_guides=wang_pos_ctrl_guides,
                                                  do_parallel=do_parallel,
                                                  draw_plots=T,
                                                  plot_dir=file.path(out_dir, "uncorrected_plots_wang"),
                                                  dep_ylims=c(-8, 3))


ggsave(file.path(out_dir, "uncorrected_HCC1419_BREAST.pdf"),
       avana_uncorrected_slopes$figs$HCC1419_BREAST,
       width=6, height=4)


p53_mutants <- ccle_mut %>%
    filter(Gene == "TP53",
           Variant_Classification != "Silent") %$%
    CellLine

p53_status <- data_frame(CellLine = unique(ccle_mut$CellLine)) %>%
    mutate(p53 = ifelse(CellLine %in% p53_mutants, "Mut", "WT"))

# p53_mutants <- ccle2maf %>%
#     filter(Hugo_Symbol == "TP53",
#            Variant_Classification != "Silent") %$%
#     Tumor_Sample_Barcode
#
# p53_status <- data_frame(CellLine = unique(ccle2maf$Tumor_Sample_Barcode)) %>%
#     mutate(p53 = ifelse(CellLine %in% p53_mutants, "Mut", "WT"))



avana_slope_plots <- slope_plots(avana_uncorrected_slopes, avana_neg_ctrl_guides, avana_pos_ctrl_guides,
                                 lineage_df, p53_status)

ggsave(file.path(out_dir, "uncorrected_fit_avana.pdf"),
       avana_slope_plots$lin_gg,
       width=3, height=3)
ggsave(file.path(out_dir, "uncorrected_fit_p53_avana.pdf"),
       avana_slope_plots$p53_gg,
       width=3, height=3)


gecko_slope_plots <- slope_plots(gecko_uncorrected_slopes, gecko_neg_ctrl_guides, gecko_pos_ctrl_guides,
                                 lineage_df, p53_status)

ggsave(file.path(out_dir, "uncorrected_fit_gecko.pdf"),
        gecko_slope_plots$lin_gg,
        width=3, height=3)
ggsave(file.path(out_dir, "uncorrected_fit_p53_gecko.pdf"),
       gecko_slope_plots$p53_gg,
       width=3, height=3)

wang_slope_plots <- slope_plots(wang_uncorrected_slopes, wang_neg_ctrl_guides, wang_pos_ctrl_guides,
                                lineage_df, p53_status, color_begin=0.4)

ggsave(file.path(out_dir, "uncorrected_fit_wang.pdf"),
       wang_slope_plots$lin_gg,
       width=3, height=3)
ggsave(file.path(out_dir, "uncorrected_fit_p53_wang.pdf"),
       wang_slope_plots$p53_gg,
       width=3, height=3)


