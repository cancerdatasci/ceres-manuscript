

slope_plots <- function(uncorrected_slopes, neg_guides, pos_guides, lineage_df,
                        p53_status=NULL,
                        x_points=seq(0,20,by=0.1),
                        min_lines=4, color_begin=0.1, color_end=0.8) {


    neg_ctrl_medians <- uncorrected_slopes$tdat %>%
        filter(Guide %in% neg_guides$Guide) %>%
        group_by(CellLine) %>%
        summarise(NegCtrl = median(Dep, na.rm=T))


    pos_ctrl_medians <- uncorrected_slopes$tdat %>%
        filter(Guide %in% pos_guides$Guide) %>%
        group_by(CellLine) %>%
        summarise(PosCtrl = median(Dep, na.rm=T))


    fit_df <- data_frame(CellLine=uncorrected_slopes$pwlm$CellLine,
                         Cuts=list(x_points)) %>% unnest(Cuts) %>%
        left_join(uncorrected_slopes$pwlm %>% select(CellLine, intercept, slope, breakpoint)) %>%
        mutate(Dep = ifelse(Cuts < breakpoint,
                            intercept + slope * Cuts,
                            intercept + slope * breakpoint)) %>%
        left_join(pos_ctrl_medians) %>%
        left_join(neg_ctrl_medians) %>%
        mutate(ScaledDep = - (Dep - NegCtrl) / (PosCtrl - NegCtrl)) %>%
        left_join(lineage_df) %>%
        group_by(Lineage) %>%
        mutate(N_Lineage = n()) %>%
        ungroup() %>%
        mutate(Lineage = ifelse(N_Lineage < 4, "Other", Lineage))


    lineage_list <- unique(fit_df$Lineage) %>% sort %>%
    {.[.!="Other"]} %>% c("Other")

    lineage_colors <- viridis::viridis(length(lineage_list),
                              begin=color_begin, end=color_end, option="B") %>%
        set_names(lineage_list)



    lin_gg <- fit_df %>% arrange(slope, Cuts) %>%
        mutate(Lineage = factor(Lineage, levels=lineage_list, ordered=T)) %>%
        ggplot(aes(x=Cuts, y=ScaledDep, group=CellLine, color=Lineage)) +
        geom_line(alpha=0.5, lwd=0.5, show.legend=F) +
        geom_hline(yintercept=c(0,-1), linetype=2, color="grey30") +
        annotate("text", size=4,
                 x=c(20, 0), y=c(0, -1), vjust=1.05, hjust=c(1.05, -0.05),
                 label=c("non-targeting controls", "core\nessentials")) +
        labs(y = "sgRNA Depletion", x="Genomic Cut Sites") +
        scale_x_continuous(limits=c(0,20), breaks=seq(0,20, by=2), expand=c(0,0)) +
        coord_cartesian(ylim=c(-1.5, 0)) +
        scale_color_manual(values=lineage_colors) +
        theme(axis.text.y=element_text(size=10),
              axis.text.x=element_text(size=9),
              axis.title.x=element_text(size=12),
              axis.title.y=element_text(size=12))




    if (!is.null(p53_status)) {
        p53_df <- fit_df %>%
            inner_join(p53_status) %>%
            group_by(p53, Cuts) %>%
            summarise(DepSD = sd(ScaledDep, na.rm=T),
                      DepMean = mean(ScaledDep, na.rm=T))

        p53_counts <- fit_df %>%
            inner_join(p53_status) %>%
            group_by(p53) %>%
            summarise(N = length(unique(CellLine)))


        p53_gg <- p53_df %>%
            left_join(p53_counts) %>%
            ungroup() %>%
            mutate(p53 = str_c(p53, " (n=", N,")")) %>%
            ggplot() +
            geom_hline(yintercept=c(0,-1), linetype=2, color="grey30") +
            geom_ribbon(aes(x=Cuts, ymin=DepMean-DepSD, ymax=DepMean+DepSD, fill=p53), alpha=0.25) +
            geom_line(aes(x=Cuts, y=DepMean, color=p53)) +
            scale_color_manual(values=c("#6a3d9a", "#33a02c")) +
            scale_fill_manual(values=c("#6a3d9a", "#33a02c")) +
            scale_x_continuous(limits=c(0,20), breaks=seq(0,20, by=2), expand=c(0,0)) +
            coord_cartesian(ylim=c(-1.5, 0)) +
            labs(y = "sgRNA Depletion", x="Genomic Cut Sites") +
            theme(legend.position=c(1, 0.95),
                  legend.justification=c(1, 1),
                  legend.title=element_text(size=12),
                  legend.text=element_text(size=12),
                  axis.text.y=element_text(size=10),
                  axis.text.x=element_text(size=9),
                  axis.title=element_text(size=12),
                  plot.title=element_text(size=12))


    } else {
        p53_gg <- NULL
    }



    return(list(lin_gg=lin_gg, p53_gg=p53_gg))


}

cn_dep_boxplot <- function(gene_level, bin_max=10, min_data_points=10,
                           ylim=c(-2, 1), boxcolors=brewer_pal(palette="Set1")(2)[1:2]) {


    cut_breaks <- c(seq(-0.5, bin_max, by=1), Inf)
    cut_labels <- 0:(bin_max-1) %>% c(str_c(bin_max,"+"))


    gene_level <- gene_level %>%
        mutate(CN_Bin = cut(CN, breaks=cut_breaks, labels=cut_labels)) %>%
        group_by(CN_Bin) %>%
        mutate(CN_Bin_N = n()) %>%
        ungroup() %>%
        filter(CN_Bin_N >= min_data_points)

    avgguide_gg <- ggplot(gene_level, aes(x=CN_Bin, y=AvgGuide)) +
        geom_hline(yintercept=c(0, -1), linetype=2, color="grey50") +
        geom_boxplot(fill=NA, outlier.size=0, outlier.stroke=0, color=boxcolors[1]) +
        scale_y_continuous(limits=ylim) +
        scale_x_discrete(drop=F) +
        labs(x = "Copy Number", y="Average Guide Score", title="Uncorrected Data")
    ceres_gg <- ggplot(gene_level, aes(x=CN_Bin, y=Gene_Effect)) +
        geom_hline(yintercept=c(0, -1), linetype=2, color="grey50") +
        geom_boxplot(fill=NA, outlier.size=0, outlier.stroke=0, color=boxcolors[2]) +
        scale_y_continuous(limits=ylim) +
        scale_x_discrete(drop=F) +
        labs(x = "Copy Number", y="CERES Gene Effect", title="CERES")

    plot_grid(avgguide_gg, ceres_gg, ncol=2)
}



cn_dep_correlation <- function(gene_level_df, title="") {
    cor_df <- gene_level_df %>% group_by(Gene) %>%
        summarise(Uncorrected = cor(AvgGuide, CN, use="pairwise"),
                  CERES = cor(Gene_Effect, CN, use="pairwise")) %>%
        gather(Method, Corr, Uncorrected, CERES)

    cor_df %>%
        ggplot(aes(x=Corr, color=Method)) +
        geom_density() +
        geom_vline(aes(xintercept=Corr, color=Method), linetype=2, show.legend=F,
                   data=cor_df %>% group_by(Method) %>%
                       summarise(Corr=mean(Corr, na.rm=T))) +
        scale_x_continuous(limits=c(-1,1), expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0)) +
        scale_color_manual(values=brewer_pal(palette="Set1")(2)[2:1]) +
        labs(x="CN-Dependency Correlation", title=title) +
        theme(axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              legend.position=c(1, 1),
              legend.justification=c(1, 1),
              legend.title=element_blank(),
              legend.text=element_text(size=10),
              plot.margin=margin(7, 14, 7, 7, "pt"))

}



breast_diff_dep_plot <- function(df, which_arms, title="",
                                 xlims=NULL, break_fdr=1e-4, max_fdr=NULL,
                                 gene_labels=NULL) {

    chr_color_scale <- scale_color_manual(values=c(chr8q="#E41A1C",
                                                   # chr8p="#377EB8",
                                                   other="black"),
                                          guide=guide_legend(override.aes=list(size=2)))

    df <- df %>%
        mutate(Qvalue = -log10(Qvalue))

    ybreaks <- seq(0, max(df$Qvalue, na.rm=T)+0.05, by=0.05)
    xbreaks <- seq(min(df$Diff)-0.1, max(df$Diff)+0.05, by=0.05)
    h2 <- hist2(df$Diff, df$Qvalue, xbreaks=xbreaks, ybreaks=ybreaks, plot=FALSE)
    xind <- cut(df$Diff, breaks=xbreaks, labels=seq(dim(h2$z)[1]), include.lowest=T, right=F)
    yind <- cut(df$Qvalue, breaks=ybreaks, labels=seq(dim(h2$z)[2]), include.lowest=T, right=F)
    df$Dens <- h2$z[cbind(xind, yind)]


    if (is.null(max_fdr)) {
        ymax <- max(df$Qvalue) + 0.05*(max(df$Qvalue) - (-log10(break_fdr)))
    } else {
        ymax <- -log10(max_fdr) + 0.05*(-log10(max_fdr) - (-log10(break_fdr)))
    }

    # if (-log10(break_fdr) > max_fdr) {
    #     break_fdr <- max_fdr
    # } else {
    #     break_fdr <- -log10(break_fdr)
    # }


    breast_chr_volcano_gg1 <-
        df %>%
        filter(!is.na(BreastCN)) %>%
        mutate(Arm = ifelse(Arm %in% which_arms, Arm, "other")) %>%
        filter(Qvalue < -log10(break_fdr)) %>%
        filter(ifelse(Arm == "other" & Dens > 1000,
                      runif(nrow(.)) < 0.005,
                      ifelse(Arm=="other" & Dens > 100,
                             runif(nrow(.)) < 0.05,
                             ifelse(Arm=="other" & Dens > 10,
                                    runif(nrow(.)) < 0.5, T)))) %>%
        mutate(GeneLabel = ifelse(Gene %in% gene_labels, Gene, "")) %>%
        # filter(ifelse(Chr == "other" & Qvalue > 0.5 & abs(Diff) < 0.05,
        #               runif(nrow(.)) < 0.001,
        #               ifelse(Chr=="other" & Qvalue < 0.5 & Qvalue > 0.25 & abs(Diff) < 0.1,
        #                      runif(nrow(.)) < 0.05, T))) %>%
        arrange(desc(Arm)) %>%
        ggplot(aes(x=Diff, y=Qvalue, color=Arm)) +
        geom_hline(yintercept=-log10(0.01), lty=2, color="grey70") +
        geom_vline(xintercept=0, lty=2, color="grey70") +
        geom_point(size=1, show.legend=F) +
        {if (!is.null(gene_labels)) {
            geom_text_repel(aes(label=GeneLabel), size=2.5,
                            show.legend=F,
                            min.segment.length=unit(3, "lines"),
                            nudge_y=0.01,
                            box.padding=unit(0.15, "lines"))
        } else {
            NULL
        }} +
        chr_color_scale +
        scale_x_continuous(limits=xlims) +
        scale_y_continuous(limits=c(0, -log10(break_fdr)), expand=c(0,0)) +
        labs(x="Differential Dependency in Breast Lines", y="-log10(FDR)")
        # theme(legend.position=c(0.95, 0.95),
        #       legend.justification=c(1, 1),
        #       legend.title=element_blank(),
        #       legend.background=element_rect(fill=NULL, size=0.5, linetype=1, color="black"))




    breast_chr_volcano_gg2 <-
        df %>%
        filter(!is.na(BreastCN)) %>%
        mutate(Arm = ifelse(Arm %in% which_arms, Arm, "other")) %>%
        filter(Qvalue >= -log10(break_fdr)) %>%
        mutate(GeneLabel = ifelse(Gene %in% gene_labels, Gene, "")) %>%
        arrange(desc(Arm)) %>%
        ggplot(aes(x=Diff, y=Qvalue, color=Arm)) +
        geom_vline(xintercept=0, lty=2, color="grey70") +
        geom_point(size=1, show.legend=T) +
        {if (!is.null(gene_labels)) {
            geom_text_repel(aes(label=GeneLabel), size=2.5,
                            show.legend=F,
                            min.segment.length=unit(3, "lines"),
                            nudge_x=0.01,
                            box.padding=unit(0.15, "lines"))
        } else {
            NULL
        }} +
        chr_color_scale +
        scale_x_continuous(limits=xlims) +
        scale_y_continuous(limits=c(-log10(break_fdr), ymax),
                                    expand=c(0,0)) +
        labs(y=" ") +
        theme(axis.line.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              legend.position=c(1, 0),
              legend.justification=c(1, 0),
              legend.title=element_blank(),
              legend.background=element_rect(fill=NULL, size=0.5, linetype=1, color="black"))

    breast_chr_dens_gg <-
        df %>%
        filter(!is.na(BreastCN)) %>%
        mutate(Arm = ifelse(Arm %in% which_arms, Arm, "other")) %>%
        mutate(Qvalue = ifelse(Qvalue < 1e-5, 1e-5, Qvalue)) %>%
        arrange(desc(Arm)) %>%
        ggplot(aes(x=Diff, color=Arm)) +
        scale_y_continuous(expand=c(0,0)) +
        chr_color_scale +
        scale_x_continuous(limits=xlims) +
        geom_density(fill=NA) +
        labs(title=title) +
        theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              plot.margin=unit(c(7,14,7,7), "pt"))

    plot_grid(breast_chr_dens_gg, breast_chr_volcano_gg2, breast_chr_volcano_gg1,
              ncol=1, align="v", rel_heights=c(1, 0.75, 3.25))
}





oncogene_plots <- function(gene_level_df, ccle_mut, gene) {

    mut_labels <- str_c(gene, c("WT", "Mut"), sep=" ")
    mut_colors <- c("black", "#ff7f00") %>% set_names(mut_labels)

    mut_df <- ccle_mut %>%
        filter(Gene == gene,
               CellLine %in% gene_level_df$CellLine)

    gene_df <- gene_level_df %>%
        filter(Gene == gene,
               CellLine %in% unique(ccle_mut$CellLine)) %>%
        left_join(mut_df) %>%
        filter(is.na(isDeleterious) | !isDeleterious) %>%
        mutate(Mut = ifelse(!is.na(Variant_Classification) &
                                Variant_Classification != "Silent",
                            mut_labels[2], mut_labels[1]) %>%
                   factor(levels=mut_labels, ordered=T),
               MutLabel = ifelse(Mut=="Mut", str_replace(Protein_Change, "p.", ""), "")) %>%
        arrange(Mut)

    uncorrected_gg <- ggplot(gene_df, aes(x=CN, y=AvgGuide, color=Mut)) +
        geom_hline(yintercept=c(0, -1), linetype=2, color="grey30") +
        geom_point() +
        labs(x="Copy Number", y=str_c("Uncorrected ", gene, " Dependency")) +
        scale_color_manual(values=mut_colors,
                           guide=guide_legend(override.aes=list(size=3))) +
        scale_x_continuous(limits=c(0, max(gene_df$CN, na.rm=T)),
                           breaks=seq(0, max(gene_df$CN, na.rm=T)+2, by=2)) +
        scale_y_continuous(breaks=seq(0, -4, by=-0.5)) +
        theme(legend.title=element_blank(),
              legend.position="right",
              legend.text=element_text(size=10),
              axis.text=element_text(size=10),
              axis.title=element_text(size=12))

    ceres_gg <- ggplot(gene_df, aes(x=CN, y=Gene_Effect, color=Mut)) +
        geom_hline(yintercept=c(0, -1), linetype=2, color="grey30") +
        geom_point() +

        labs(x="Copy Number", y=str_c("CERES ", gene, " Dependency")) +
        scale_color_manual(values=mut_colors,
                           guide=guide_legend(override.aes=list(size=3))) +
        scale_x_continuous(limits=c(0, max(gene_df$CN, na.rm=T)),
                           breaks=seq(0, max(gene_df$CN, na.rm=T)+2, by=2)) +
        scale_y_continuous(breaks=seq(0, -4, by=-0.5)) +
        theme(legend.title=element_blank(),
              legend.position="right",
              legend.text=element_text(size=10),
              axis.text=element_text(size=10),
              axis.title=element_text(size=12))

    return(list(uncorrected=uncorrected_gg, ceres=ceres_gg))

}
