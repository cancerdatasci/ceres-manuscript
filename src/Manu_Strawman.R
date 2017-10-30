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

out_dir <- file.path("./output/figures/strawman_correction", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")


load("./cache/ne_genes.RData")
load("./cache/cce_genes.RData")


dep_cn_lm <- avana_gene_level %>%
    select(Gene, CellLine, CN, Dep=AvgGuide) %>%
    group_by(CellLine) %>%
    do(broom::augment(lm(.$Dep ~ .$CN), newdata=.)) %>%
    mutate(Residuals = Dep - .fitted)


avana_strawman_genelevel <-
    make.a.tidy.dataset(list(LM_Dep = dep_cn_lm %>%
                                 select(Gene, CellLine, Residuals) %>%
                                 df.to.mat(),
                             CERES = avana_gene_level %>%
                                 select(Gene, CellLine, Gene_Effect) %>%
                                 df.to.mat,
                             Uncorrected = avana_gene_level %>%
                                 select(Gene, CellLine, AvgGuide) %>%
                                 df.to.mat),
                        dim.names=c("Gene", "CellLine"))


avana_raw_lm_genelevel_pr <-
    run_precision_recall(avana_strawman_genelevel, neg_genes=ne_genes, pos_genes=cce_genes,
                         case0_col="Uncorrected", case0_label="Uncorrected Data",
                         case1_col="LM_Dep", case1_label="LM Correction", x_axis_pos="bottom")


avana_raw_lm_genelevel_pr$gg +
    theme(legend.position=c(1,0), legend.justification=c(1,0)) +
    ggsave(file.path(out_dir, "strawman_lm_p-r_avana.pdf"),
           width=8, height=4)


cn_thresh <- 4
dep_thresh <- -0.6
filtering_amps <-
    avana_gene_level %>%
    filter(CN > cn_thresh) %>%
    group_by(CellLine) %>%
    summarise(TotalAmps = n(),
              EssentialAmps = sum(Gene_Effect < dep_thresh))

ggplot(filtering_amps, aes(x=TotalAmps)) +
    geom_histogram(binwidth=20, fill="grey30") +
    geom_vline(xintercept = mean(filtering_amps$TotalAmps), linetype=2, color="grey30") +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(x = "Number of Genes Filtered", y = "Number of Cell Lines",
         title = "All Genes") +
    ggsave(file.path(out_dir, "strawman_filtered_all_avana.pdf"),
           width=4, height=3)

ggplot(filtering_amps, aes(x=EssentialAmps)) +
    geom_histogram(binwidth=1, fill="#E41A1C") +
    geom_vline(xintercept = mean(filtering_amps$EssentialAmps), linetype=2, color="#E41A1C") +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    labs(x = "Number of Genes Filtered", y = "Number of Cell Lines",
         title = "Essential Genes") +
    ggsave(file.path(out_dir, "strawman_filtered_essential_avana.pdf"),
           width=4, height=3)


