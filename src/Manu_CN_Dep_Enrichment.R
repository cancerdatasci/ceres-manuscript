library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/cn_dep_enrichment", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


load("./cache/lineage_df.RData")

avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")
gecko_gene_level <- readRDS("./data/ceres_cache_final/gecko_gene_level.rds")
wang_gene_level <- readRDS("./data/ceres_cache_final/wang_gene_level.rds")


top_n_genes <- 100

avana_amp_enrichment <-
    run_cn_dep_enrichment(avana_gene_level, lineage_df, top_n_genes,
                                 dep_col=AvgGuide, cn_col=CN)


ggsave(file.path(out_dir, "dep_cn_top100_avana.pdf"),
       avana_amp_enrichment$gg,
       width=10, height=4)

ggsave(file.path(out_dir, "dep_cn_top100_ks_avana.pdf"),
       avana_amp_enrichment$ks_gg,
       width=4, height=4)


gecko_amp_enrichment <-
    run_cn_dep_enrichment(gecko_gene_level, lineage_df, top_n_genes,
                          dep_col=AvgGuide, cn_col=CN)


ggsave(file.path(out_dir, "dep_cn_top100_gecko.pdf"),
       gecko_amp_enrichment$gg,
       width=8, height=4)

ggsave(file.path(out_dir, "dep_cn_top100_ks_gecko.pdf"),
       gecko_amp_enrichment$ks_gg,
       width=4, height=4)

wang_amp_enrichment <-
    run_cn_dep_enrichment(wang_gene_level, lineage_df, top_n_genes,
                          dep_col=AvgGuide, cn_col=CN)

ggsave(file.path(out_dir, "dep_cn_top100_wang.pdf"),
       wang_amp_enrichment$gg,
       width=5, height=4)

ggsave(file.path(out_dir, "dep_cn_top100_ks_wang.pdf"),
       wang_amp_enrichment$ks_gg,
       width=4, height=4)




cn_breaks <- c(4, 6, 8, Inf)
cn_labels <- c("4-6", "6-8", "8+") %>% str_c("CN: ", .)
cn_colors <- brewer_pal(palette="YlOrRd")(9)[c(5,7,9)]


essentials <- c(
    "KEGG_DNA_REPLICATION" = "DNA Replication",
    "KEGG_RNA_POLYMERASE" = "RNA Polymerase",
    "KEGG_RIBOSOME" = "Ribosome",
    "KEGG_SPLICEOSOME" = "Spliceosome",
    "KEGG_PROTEASOME" = "Proteasome")

essential_gene_sets <- read.gmt("./data/raw/c2.all.v6.0.symbols.gmt") %>%
    filter(GeneSet %in% names(essentials)) %>%
    mutate(EssentialSet = essentials[as.character(GeneSet)])

fp_df <- avana_gene_level %>% group_by(CellLine) %>%
    mutate(DepRank = min_rank(AvgGuide),
           CN_Bin = cut(CN, cn_breaks, cn_labels)) %>%
    left_join(lineage_df) %>%
    group_by(Lineage) %>%
    mutate(LineageN = length(unique(CellLine))) %>%
    ungroup() %>%
    mutate(Lineage = ifelse(LineageN < 4, "Other", Lineage))

lineage_list <- unique(fp_df$Lineage) %>% sort %>%
{.[.!="Other"]} %>% c("Other")

lineage_colors <- viridis::viridis(length(lineage_list),
                                   begin=0.1, end=0.8, option="B") %>%
    set_names(lineage_list)

example_cell_line <- "HCC1419_BREAST"

cell_line_color <- lineage_colors[lineage_df %>%
                                      filter(CellLine==example_cell_line) %$%
                                      Lineage]


x_lims <- range(fp_df$DepRank, na.rm=T)

cell_line_label <- str_extract(example_cell_line, "^[^_]+")

cl_fp_df <- fp_df %>%
    filter(CellLine == example_cell_line)

waterfall_gg <- cl_fp_df %>%
    mutate(Dummy = "Average Guide Score") %>%
    ggplot(aes(x=DepRank, y=AvgGuide)) +
    geom_hline(yintercept=c(0, -1), linetype=2, color="grey30") +
    facet_grid( Dummy ~ ., switch="y")  +
    geom_line(aes(color=Lineage), size=1.5, show.legend=F) +
    scale_x_continuous(limits=x_lims,
                       breaks=c(1,5000,10000,15000),
                       expand=c(0.01, 0),
                       position="top") +
    scale_y_continuous(breaks=pretty_breaks()) +
    scale_color_manual(values=lineage_colors) +
    labs(x = "Genes Ranked by Depletion", y = "Average Guide Score") +
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          strip.placement="outside",
          axis.title.y=element_blank())

geneset_df <- bind_rows(cl_fp_df %>%
                            filter(Gene %in% essential_gene_sets$Gene) %>%
                            left_join(essential_gene_sets) %>%
                            select(CellLine, Gene, DepRank, Set = EssentialSet),
                        cl_fp_df %>%
                            filter(!is.na(CN_Bin)) %>%
                            select(CellLine, Gene, DepRank, Set = CN_Bin)) %>%
    mutate(Set = factor(Set, levels=c(essentials, cn_labels), ordered=T))

geneset_gg <-
    ggplot(geneset_df) +
    facet_grid(Set ~ ., scales="free", switch="y") +
    geom_point(aes(y=Set, x=-1)) +
    geom_vline(aes(xintercept=DepRank, color=Set),
               alpha=0.5, size=0.5, show.legend=F) +
    scale_color_manual(values=c(cn_colors %>% set_names(cn_labels),
                                rep("black",length(essentials)) %>% set_names(essentials))) +
    scale_x_continuous(limits=x_lims, expand=c(0, 0)) +
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          strip.placement="outside",
          panel.spacing=unit(0,"pt"),
          panel.background=element_rect(linetype=1, size=0.5, color="grey20"),
          panel.ontop=T,
          axis.line=element_blank(),
          axis.ticks=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=12),
          axis.title=element_blank())

plot_grid(waterfall_gg, geneset_gg, ncol=1, align="v",
          rel_heights=c(1, 1)) +
    ggsave(file.path(out_dir, str_c(example_cell_line, ".pdf")),
           width=4, height=4)


