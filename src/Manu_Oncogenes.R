library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/oncogenes", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")


load("./cache/lineage_df.RData")
load("./cache/ccle_mut.RData")

kras_avana_plots <- oncogene_plots(avana_gene_level, ccle_mut, "KRAS")

ggsave(file.path(out_dir, "kras_uncorrected_avana.pdf"),
       kras_avana_plots$uncorrected, width=4, height=3)

ggsave(file.path(out_dir, "kras_ceres_avana.pdf"),
       kras_avana_plots$ceres, width=4, height=3)

braf_avana_plots <- oncogene_plots(avana_gene_level, ccle_mut, "BRAF")

ggsave(file.path(out_dir, "braf_uncorrected_avana.pdf"),
       braf_avana_plots$uncorrected, width=4, height=3)

ggsave(file.path(out_dir, "braf_ceres_avana.pdf"),
       braf_avana_plots$ceres, width=4, height=3)


pik3ca_avana_plots <- oncogene_plots(avana_gene_level, ccle_mut, "PIK3CA")

ggsave(file.path(out_dir, "pik3ca_uncorrected_avana.pdf"),
       pik3ca_avana_plots$uncorrected, width=4, height=3)

ggsave(file.path(out_dir, "pik3ca_ceres_avana.pdf"),
       pik3ca_avana_plots$ceres, width=4, height=3)



lin_labels <- c("Neuroblastoma", "Other")
lin_colors <- c("#984ea3", "black") %>% set_names(lin_labels)

nb_df <- avana_gene_level %>%
    filter(Gene == "MYCN") %>%
    left_join(lineage_df) %>%
    mutate(Lineage = ifelse(Lineage=="Neuroblastoma",
                            "Neuroblastoma", "Other") %>%
               factor(levels=lin_labels, ordered=T)) %>%
    arrange(desc(Lineage))


mycn_uncorrected_gg <- ggplot(nb_df, aes(x=CN, y=AvgGuide, color=Lineage)) +
    geom_hline(yintercept=c(0, -1), linetype=2, color="grey30") +
    geom_point() +
    labs(x="Copy Number", y=str_c("Uncorrected MYCN Dependency")) +
    scale_color_manual(values=lin_colors,
                       guide=guide_legend(override.aes=list(size=3))) +
    scale_x_continuous(limits=c(0, max(nb_df$CN, na.rm=T)),
                       breaks=seq(0, max(nb_df$CN, na.rm=T)+4, by=4)) +
    scale_y_continuous(breaks=seq(0, -4, by=-0.5)) +
    theme(legend.title=element_blank(),
          legend.position="top",
          legend.text=element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))

mycn_ceres_gg <- ggplot(nb_df, aes(x=CN, y=Gene_Effect, color=Lineage)) +
    geom_hline(yintercept=c(0, -1), linetype=2, color="grey30") +
    geom_point() +

    labs(x="Copy Number", y=str_c("CERES MYCN Dependency")) +
    scale_color_manual(values=lin_colors,
                       guide=guide_legend(override.aes=list(size=3))) +
    scale_x_continuous(limits=c(0, max(nb_df$CN, na.rm=T)),
                       breaks=seq(0, max(nb_df$CN, na.rm=T)+4, by=4)) +
    scale_y_continuous(breaks=seq(0, -4, by=-0.5)) +
    theme(legend.title=element_blank(),

          legend.position="top",
          legend.text=element_text(size=10),
          axis.text=element_text(size=10),
          axis.title=element_text(size=12))


ggsave(file.path(out_dir, "mycn_uncorrected_avana.pdf"),
       mycn_uncorrected_gg, width=4, height=3)

ggsave(file.path(out_dir, "mycn_ceres_avana.pdf"),
       mycn_ceres_gg, width=4, height=3)
