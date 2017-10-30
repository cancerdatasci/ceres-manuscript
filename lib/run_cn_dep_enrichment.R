
run_cn_dep_enrichment <- function(df, lineage_df, top_n_genes=100,
                                  dep_col=AvgGuide, cn_col=CN,
                                  cn_breaks=c(2, 4, 6, 8, Inf),
                                  cn_labels=c("2-4  ", "4-6  ", "6-8  ", "8+  "),
                                  cn_colors=brewer_pal(palette="YlOrRd")(9)[c(4,5,7,9)],
                                  min_lineage_group=4) {

    dep_col <- enquo(dep_col)
    cn_col <- enquo(cn_col)



    rank_df <- df %>%
        filter(!is.na(!!cn_col), !is.na(!!dep_col)) %>%
        group_by(CellLine) %>%
        mutate(Rank = min_rank(!!dep_col),
               TopN = min_rank(-!!cn_col) <= top_n_genes,
               Total = length(unique(Gene))) %>%
        mutate(ScaledRank = Rank/Total)


    cell_line_top_amps <-
        rank_df %>% group_by(CellLine) %>%
        summarise(MeanCN = mean((!!cn_col)[TopN]),
                  MedianRank = median(ScaledRank[TopN]),
                  Quart1Rank = quantile(ScaledRank[TopN], 0.25, na.rm=T),
                  Quart3Rank = quantile(ScaledRank[TopN], 0.75, na.rm=T)) %>%
        ungroup() %>%
        left_join(lineage_df) %>%
        group_by(Lineage) %>%
        mutate(LineageN = n()) %>%
        ungroup() %>%
        mutate(Lineage = ifelse(LineageN < min_lineage_group, "Other", Lineage))

    lineage_list <- sort(unique(cell_line_top_amps$Lineage)) %>%
    {.[.!="Other"]} %>% c("Other")


    # ylims <- range()
    # ylims <- c(1, df %>% filter(!is.na(CN), !is.na(Dep)) %$% length(unique(Gene)))
    ylims <- c(0,1)

    # cn_breaks <- c(2, 4, 6, 8, Inf)
    # cn_labels <- c("2-4  ", "4-6  ", "6-8  ", "8+  ")
    # cn_colors <- brewer_pal(palette="YlOrRd")(9)[c(4,5,7,9)]


    tmp <- cell_line_top_amps %>%
        mutate(CN = cut(MeanCN, breaks=cn_breaks, labels=cn_labels)) %>%
        mutate(Lineage = factor(Lineage, levels=lineage_list, ordered=T)) %>%
        arrange(Lineage, MedianRank) %>%
        # group_by(Lineage) %>%
        mutate(CellLine = row_number() + (as.numeric(Lineage)-1)*10)




    lineage_labs <- tmp %>%
        group_by(Lineage) %>%
        summarise(CellLineMin = min(CellLine),
                  CellLineMid = mean(CellLine)) %$%

        #rep("", length(CellLine)) %>%
        #              inset(floor(CellLineMid-CellLineMin+1), as.character(Lineage[1]))) %>%
        #ungroup() filter(CellLineLabs != "") %$%
        set_names(Lineage, CellLineMid)

    gg <- tmp %>%
        # bind_rows(.,
        #           distinct(., Lineage, CellLineMin) %>%
        #               mutate(CellLine=CellLineMin-0.001, MedianRank=-1)) %>%

        ggplot(aes(x=CellLine, y=MedianRank, color=CN)) +
        facet_grid( ~ Lineage, scales="free_x", space="free_x", switch="x") +
        geom_rect(aes(xmin=-Inf, xmax=Inf,
                      ymin=ylims[2]*0.75, ymax=ylims[2]*0.25),
                  fill="grey95", inherit.aes=F) +
        geom_hline(yintercept = ylims[2]*c(0.25, 0.5, 0.75), color="grey20", size=0.25) +
        geom_segment(aes(x=CellLine, xend=CellLine, y=Quart1Rank, yend=Quart3Rank),
                     size=0.3) +
        geom_point(size=0.5, show.legend=F) +
        scale_x_continuous(expand=c(0, 1.5),
                           position="top",
                           breaks=names(lineage_labs) %>% as.numeric(),
                           labels=lineage_labs) +
        scale_y_continuous(limits=ylims,
                           # breaks=c(1, 5000, 10000, 15000),
                           breaks=c(0, 0.25, 0.5, 0.75, 1),
                           expand=c(0, 0)) +
        # scale_color_continuous(trans="log10", name="Number of amplified genes") +
        scale_color_manual(values=cn_colors %>% set_names(cn_labels),
                           name="Mean Copy Number:", na.translate=F,
                           guide=guide_legend(override.aes=list(size=c(1.5)))) +
        labs(y = "Depletion Scaled Rank\n(Interquartile Range)") +
        theme(plot.title=element_text(size=12, face="plain"),
              legend.title=element_text(size=12),
              legend.text=element_text(size=12),
              legend.position="bottom",
              legend.margin=margin(0, 0, 0, 0, "pt"),
              axis.title.y=element_text(size=12),
              axis.text.y=element_text(size=8),
              axis.text.x=element_text(size=10, angle=-30, vjust=1, hjust=1),
              axis.title.x=element_blank(),
              axis.ticks.x=element_line(size=0.25),
              axis.ticks.y=element_line(size=0.25),
              axis.line.x=element_blank(),
              axis.line.y=element_blank(),
              strip.background=element_blank(),
              # strip.text=element_text(size=8, angle=90, vjust=1),
              strip.text=element_blank(),
              panel.ontop=T,
              panel.spacing=unit(2, "pt"),
              panel.background=element_rect(size=0.5, linetype=1,
                                            color="grey20", fill=NULL))

    #
    # ceres_ranks <- avana_gene_level %>%
    #     filter(!is.na(CN), !is.na(AvgGuide), !is.na(Gene_Effect)) %>%
    #     distinct(CellLine, AvgGuide, .keep_all=T) %>%
    #     distinct(CellLine, Gene_Effect, .keep_all=T) %>%
    #     group_by(CellLine) %>%
    #     mutate(Rank = min_rank(Gene_Effect),
    #            Total = length(unique(Gene)),
    #            TopN = min_rank(-CN) <= top_n) %>%
    #     filter(TopN)



    punif_discrete <- function(q, min, max) {(q - min + 1) / (max - min + 1)}

    ks_df <- rank_df %>%
        filter(TopN) %>%
        group_by(CellLine) %>%
        mutate(MeanCN = mean(CN)) %>%
        group_by(CellLine, MeanCN) %>%
        do(broom::tidy(ks.test(.$Rank, "punif_discrete", min=1, max=.$Total[1],
                               alternative="greater"))) %>%
        ungroup %>%
        mutate(q.value = p.adjust(p.value, method="fdr"))

    ks_gg <- ggplot(ks_df, aes(x=MeanCN, y=statistic, color=p.value < 0.05)) +
        geom_point(size=1, show.legend=F) +
        labs(x = "Mean CN", y = "K-S Enrichment Statistic",
             title="100 Most Amplified Genes") +
        scale_color_manual(values=c("FALSE"="black", "TRUE"="#e41a1c"))

    return(list(gg=gg, ks_df=ks_df, ks_gg=ks_gg))
}
