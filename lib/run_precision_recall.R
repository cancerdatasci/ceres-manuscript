
run_precision_recall <- function(gene_level_df, pos_genes, neg_genes, fdr=5,
                                 xlab="Cell Line", x_axis_pos="top",
                                 case0_col="AvgGuide", case1_col="Gene_Effect",
                                 case0_label="Uncorrected Data", case1_label="CERES") {

    gene_level_df <-
        gene_level_df %>%
        dplyr::rename_(Case0=case0_col,
                       Case1=case1_col)

    case0_precision_recall <- pr_curve(gene_level_df, Case0, neg_genes, pos_genes)
    # case0_precision_recall <- gene_level_df %>%
    #     filter(Gene %in% c(neg_genes, pos_genes)) %>%
    #     mutate(Group = ifelse(Gene %in% neg_genes, "NE", "CCE")) %>%
    #     arrange(CellLine, Case0) %>%
    #     group_by(CellLine) %>%
    #     mutate(CCE_Sum = cumsum(Group == "CCE"),
    #            NE_Sum = cumsum(Group == "NE")) %>%
    #     mutate(Recall = CCE_Sum / sum(Group == "CCE"),
    #            Precision = CCE_Sum / (CCE_Sum + NE_Sum)) %>%
    #     select(CellLine, Gene, Recall, Precision, Group)

    case1_precision_recall <- pr_curve(gene_level_df, Case1, neg_genes, pos_genes)
    # case1_precision_recall <- gene_level_df %>%
    #     filter(Gene %in% c(neg_genes, pos_genes)) %>%
    #     mutate(Group = ifelse(Gene %in% neg_genes, "NE", "CCE")) %>%
    #     arrange(CellLine, Case1) %>%
    #     group_by(CellLine) %>%
    #     mutate(CCE_Sum = cumsum(Group == "CCE"),
    #            NE_Sum = cumsum(Group == "NE")) %>%
    #     mutate(Recall = CCE_Sum / sum(Group == "CCE"),
    #            Precision = CCE_Sum / (CCE_Sum + NE_Sum)) %>%
    #     select(CellLine, Gene, Recall, Precision, Group)

    case0_fdr_recall <- recall_at_fdr(case0_precision_recall, fdr)
    # case0_fdr_recall <- case0_precision_recall %>%
    #     group_by(CellLine) %>%
    #     summarise(Recall = max(c(0, Recall[Precision >= 1-fdr/100])))

    case1_fdr_recall <- recall_at_fdr(case1_precision_recall, fdr)
    # case1_fdr_recall <- case1_precision_recall %>%
    #     group_by(CellLine) %>%
    #     summarise(Recall = max(c(0, Recall[Precision >= 1-fdr/100])))

    recall_df <- inner_join(case0_fdr_recall %>% dplyr::rename(Case0_Recall = Recall),
                            case1_fdr_recall %>% dplyr::rename(Case1_Recall = Recall)) %>%
        mutate(Avg = (Case0_Recall + Case1_Recall) / 2) %>%
        # arrange(Avg) %>%
        arrange(Case0_Recall, Case1_Recall) %>%
        mutate(CellLine = factor(CellLine, levels=CellLine, ordered=T))

    recall_gather_df <- recall_df %>%
        gather(Method, Recall, -CellLine, -Avg) %>%
        mutate(Method = factor(Method,
                               levels=c("Case1_Recall", "Case0_Recall"), ordered=T))
    # %>%
    #     mutate(Method = factor(Method, levels=c("Uncorrected Data", "CERES"),
    #                            ordered=T))

    gg <- ggplot(recall_gather_df) +
        geom_segment(aes(x=CellLine, xend=CellLine, y=Case0_Recall, yend=Case1_Recall),
                     data=recall_df, color="grey80", size=0.5) +
        geom_point(aes(x=CellLine, y=Recall, color=Method)) +
        scale_y_continuous(limits=c(-0.01,1), breaks=seq(0, 1, 0.2), expand=c(0,0)) +
        scale_x_discrete(expand=c(0, 2), position=x_axis_pos) +
        scale_color_manual(values=brewer_pal(palette="Set1")(2)[2:1],
                           labels=c(Case0_Recall=case0_label, Case1_Recall=case1_label),
                           guide=guide_legend(override.aes=list(size=3))) +
        labs(x=xlab, y=str_c("Recall at ", fdr, "% FDR")) +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              # legend.position="top",
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.position=c(0, 1),
              legend.justification=c(0, 1))

    return(list(df=recall_df, gg=gg, case0=case0_precision_recall, case1=case1_precision_recall))
}


pr_curve <- function(df, col, neg_genes, pos_genes) {

    col <- enquo(col)

    pr <- df %>%
        filter(Gene %in% c(neg_genes, pos_genes)) %>%
        mutate(Group = ifelse(Gene %in% neg_genes, "NE", "CCE")) %>%
        arrange(CellLine, !!col) %>%
        group_by(CellLine) %>%
        mutate(CCE_Sum = cumsum(Group == "CCE"),
               NE_Sum = cumsum(Group == "NE")) %>%
        mutate(Recall = CCE_Sum / sum(Group == "CCE"),
               Precision = CCE_Sum / (CCE_Sum + NE_Sum)) %>%
        select(CellLine, Gene, Recall, Precision, Group)
}

recall_at_fdr <- function(pr, fdr) {
    recall <- pr %>%
        group_by(CellLine) %>%
        summarise(Recall = max(c(0, Recall[Precision >= 1-fdr/100])))
}

