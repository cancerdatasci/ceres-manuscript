
run_auc <- function(gene_level, lineage_df, ne_genes, cce_genes) {

    roc_df <- gene_level %>%
        filter(!is.na(CN), !is.na(AvgGuide), !is.na(Gene_Effect)) %>%
        mutate(CCE = Gene %in% cce_genes,
               NE = Gene %in% ne_genes) %>%
        filter(CCE | NE) %>%
        arrange(CellLine, AvgGuide) %>%
        group_by(CellLine) %>%
        mutate(TPR = cumsum(CCE) / sum(CCE),
               FPR = cumsum(NE) / sum(NE))

    auc_df <- roc_df %>%
        group_by(CellLine) %>%
        summarise(AUC = sum(TPR[-1] * diff(FPR))) %>%
        ungroup %>%
        left_join(lineage_df) %>%
        group_by(Lineage) %>%
        mutate(N_Lineage = n()) %>%
        ungroup
}



auc_sunburst <- function(auc_df, min_lines=4, min_label_size=2.5, max_label_size=4.5,
                         label_pos_y=0.89, label_pos_x=0.5,
                         fill_begin=0.1, fill_end=0.8) {

    sunburst_df <- auc_df %>%
        mutate(Lineage = ifelse(N_Lineage < min_lines, "Other", Lineage),
               CellLine = str_extract(CellLine, "^[^_]+")) %>%
        mutate(Lineage = factor(Lineage, levels=sort(unique(auc_df$Lineage)) %>% c("Other"),
                                ordered=T)) %>%
        arrange(Lineage, -AUC) %>% ungroup %>%
        mutate(CellLine = factor(CellLine, levels=CellLine, ordered=T),
               TextAngle = ifelse(1:n() < n()/2,
                                  seq(90, -270, length.out = n()),
                                  180 + seq(90, -270, length.out = n())),
               HJust = ifelse(1:n() < n()/2, 1, 0))


    sunburst_lineages <- sunburst_df %>%
        group_by(Lineage) %>%
        summarise(X = label_pos_x*(max(as.numeric(CellLine))-min(as.numeric(CellLine))) +
                      min(as.numeric(CellLine)),
                  N = n()) %>%
        ungroup() %>%
        mutate(TextAngle = rescale(X, from=c(0.5,nrow(sunburst_df)+0.5), to=c(90, -270)) +
                   ifelse(X < nrow(sunburst_df)/2, 0, 180),
               HJust = ifelse(X < nrow(sunburst_df)/2, 0, 1))

    lineage_list <- unique(sunburst_lineages$Lineage) %>% sort

    lineage_colors <- viridis::viridis(length(lineage_list),
                              begin=fill_begin, end=fill_end, option="B") %>%
        set_names(lineage_list)


    bind_rows(sunburst_df %>% mutate(CellLine = as.numeric(CellLine)-0.5),
              sunburst_df %>% mutate(CellLine = as.numeric(CellLine)+0.5)) %>%
        arrange(Lineage, CellLine, -AUC) %>%
        ggplot() +

        geom_ribbon(aes(x=CellLine, ymin=0, ymax=AUC-0.5, fill=Lineage), show.legend=F, color="white") +
        geom_hline(yintercept=c(0.3, 0.4, 0.5), color="black", linetype=c(3, 3, 1)) +
        geom_text(aes(x=X, label=str_c(Lineage, " (", N, ")"),
                      angle=TextAngle, hjust=HJust, y=label_pos_y-0.5,
                      size=sqrt(N)), vjust=0.5, show.legend=F,
                  color="white", fontface="bold", data=sunburst_lineages) +
        annotate(geom="text", x=(nrow(sunburst_df)+1)/2, y=c(0.3, 0.4, 0.5), label=c("0.8", "0.9", "AUC=1.0"),
                 vjust=1.25, colour="black", size=6, fontface="bold") +
        coord_polar(start=pi) +
        scale_x_continuous(limits=c(0.5, nrow(sunburst_df)+0.5), expand=c(0,0)) +
        scale_fill_manual(values=lineage_colors) +
        scale_size_continuous(range=c(min_label_size, max_label_size)) +
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              plot.margin=margin(-4, -4, -4, -4, "lines"))

}
