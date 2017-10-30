library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/precision_recall", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)

load("./cache/ne_genes.RData")
load("./cache/cce_genes.RData")
avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")
gecko_gene_level <- readRDS("./data/ceres_cache_final/gecko_gene_level.rds")
wang_gene_level <- readRDS("./data/ceres_cache_final/wang_gene_level.rds")


avana_precision_recall <-
    run_precision_recall(avana_gene_level, cce_genes, ne_genes, x_axis_pos="top",
                         case0_label="Uncorrected Data", case1_label="CERES")
avana_precision_recall$gg +
    theme(legend.direction="horizontal") +
    ggsave(file.path(out_dir, "avana_cellline_recall.pdf"), width=8, height=6)

gecko_precision_recall <-
    run_precision_recall(gecko_gene_level, cce_genes, ne_genes, x_axis_pos="bottom",
                         case0_label="Uncorrected Data", case1_label="CERES")
gecko_precision_recall$gg +
    theme(legend.position=c(1,0),
          legend.justification=c(1,0)) +
    ggsave(file.path(out_dir, "gecko_cellline_recall.pdf"), width=4, height=4)

wang_precision_recall <-
    run_precision_recall(wang_gene_level, cce_genes, ne_genes, x_axis_pos="bottom",
                         case0_label="Uncorrected Data", case1_label="CERES")
wang_precision_recall$gg +
    theme(legend.position=c(1,0),
          legend.justification=c(1,0)) +
    ggsave(file.path(out_dir, "wang_cellline_recall.pdf"), width=4, height=4)



p_r_dir <- file.path(out_dir, "precision_recall_curves")
dir.create(p_r_dir, recursive=T, showWarnings=F)



precision_recall_curves <- function(gene_level_df, pos_genes, neg_genes, fdr=5) {

    avgguide_precision_recall <- gene_level_df %>%
        filter(Gene %in% c(neg_genes, pos_genes)) %>%
        mutate(Group = ifelse(Gene %in% neg_genes, "NE", "CCE")) %>%
        arrange(CellLine, AvgGuide) %>%
        group_by(CellLine) %>%
        mutate(CCE_Sum = cumsum(Group == "CCE"),
               NE_Sum = cumsum(Group == "NE")) %>%
        mutate(Recall = CCE_Sum / sum(Group == "CCE"),
               Precision = CCE_Sum / (CCE_Sum + NE_Sum)) %>%
        select(CellLine, Gene, Dep=AvgGuide, Recall, Precision, Group) %>%
        mutate(Dataset = "Uncorrected Data")

    ceres_precision_recall <- gene_level_df %>%
        filter(Gene %in% c(neg_genes, pos_genes)) %>%
        mutate(Group = ifelse(Gene %in% neg_genes, "NE", "CCE")) %>%
        arrange(CellLine, Gene_Effect) %>%
        group_by(CellLine) %>%
        mutate(CCE_Sum = cumsum(Group == "CCE"),
               NE_Sum = cumsum(Group == "NE")) %>%
        mutate(Recall = CCE_Sum / sum(Group == "CCE"),
               Precision = CCE_Sum / (CCE_Sum + NE_Sum)) %>%
        select(CellLine, Gene, Dep=Gene_Effect, Recall, Precision, Group) %>%
        mutate(Dataset = "CERES")


    avgguide_fdr_recall <- avgguide_precision_recall %>%
        group_by(CellLine) %>%
        summarise(UncorrectedRecall = max(c(0, Recall[Precision >= 1-fdr/100])))


    ceres_fdr_recall <- ceres_precision_recall %>%
        group_by(CellLine) %>%
        summarise(CERESRecall = max(c(0, Recall[Precision >= 1-fdr/100])))

    fdr_recall <- inner_join(avgguide_fdr_recall, ceres_fdr_recall) %>%
        arrange(UncorrectedRecall, CERESRecall) %>%
        mutate(Bin = row_number() %>% cut(breaks=seq(0, n()+5, by=5))) %>%
        mutate(Bin = factor(Bin) %>% as.numeric)




    ggs <- bind_rows(avgguide_precision_recall, ceres_precision_recall) %>%
        left_join(fdr_recall) %>%
        dlply(.(CellLine), function(df) {


            cell_line <- df$CellLine[1]
            cell_line_short <- str_extract(cell_line, "^[^_]+")



            df %>%
                ggplot(aes(x=Recall, y=Precision, color=Dataset)) +
                geom_hline(yintercept = 1-fdr/100,
                           color="grey30", linetype=2) +
                geom_line(show.legend=F) +
                annotate("text", x=0, y=1-fdr/100, hjust=-0.1, vjust=-0.2,
                         label=str_c(fdr, "% FDR"), size=3) +
                scale_x_continuous(limits=c(0,1), expand=c(0,0),
                                   breaks=seq(0, 1, 0.2)) +
                coord_cartesian(ylim=c(0.7, 1)) +
                scale_color_manual(values=brewer_pal(palette="Set1")(2)[2:1]) +
                labs(title=cell_line_short) +
                theme(legend.title=element_blank(),
                      legend.margin=margin(0,0,0,0, "pt"),
                      legend.text=element_text(size=8),
                      axis.text=element_text(size=8),
                      axis.title=element_text(size=10),
                      plot.title=element_text(size=10)) +
                ggsave(file.path(p_r_dir, str_c(cell_line, ".pdf")),
                       width=2, height=2)
        })
}

# cls_to_plot <-
#     c(avana_precision_recall$df %>%
#           top_n(50, Case0_Recall) %$% as.character(CellLine),
#       avana_precision_recall$df %>%
#           top_n(50, -Case0_Recall) %$% as.character(CellLine))

ggs <- avana_gene_level %>%
    # filter(CellLine %in% cls_to_plot) %>%
    precision_recall_curves(cce_genes, ne_genes, 5)



