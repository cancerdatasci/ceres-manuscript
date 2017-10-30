library(ggrepel)

library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))

threads <- 2

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- T
} else {
    do_parallel <- F
}

out_dir <- file.path("./output/figures/lineage_arm_deps", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)

load("./cache/avana_genes.RData")
load("./cache/xpr.RData")
load("./cache/cn_genes.RData")
load("./cache/lineage_df.RData")

avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")


rnai_dep <-
    make.a.tidy.dataset(list(RNAi = read.gct("./data/raw/RNAi_demeter_gene_deps.gct") %>%
                                 set_rownames(str_extract(rownames(.), "^[^\\s]+")),
                             RPKM = xpr %>% select(Gene, CellLine, RPKM) %>%
                                 distinct(Gene, CellLine, .keep_all=T) %>% df.to.mat(),
                             CN = cn_genes %>% select(Gene, CellLine, CN) %>% df.to.mat()),
                        use.dims="RNAi", dim.names=c("Gene", "CellLine")) %>%
    left_join(lineage_df)


cytobands <- read_tsv("./data/raw/hg19_cytoBand.txt",
                      col_names=c("Chr", "Start", "End", "Band", "Stain"))

arm_coords <- cytobands %>%
    mutate(ArmLabel = str_sub(Band, 1, 1)) %>%
    group_by(Chr, ArmLabel) %>%
    summarise(Start = min(Start) + 1,
              End = max(End)) %>%
    ungroup %>%
    mutate(Arm = str_c(Chr, ArmLabel))


gene_arms <-
    avana_genes %>% distinct(Gene, Chr, CDS_Start, CDS_End) %>%
    left_join(arm_coords) %>%
    filter(CDS_Start < End, CDS_End > Start)


breast_v_rest_avgguide <-
    avana_gene_level %>%
    mutate(Breast = str_detect(CellLine, "_BREAST")) %>%
    group_by(Gene) %>%
    summarise(BreastMean = mean(AvgGuide[Breast], na.rm=T),
              OtherMean = mean(AvgGuide[!Breast], na.rm=T),
              Diff = BreastMean - OtherMean,
              Pvalue = t.test(AvgGuide[Breast], AvgGuide[!Breast],
                              var.equal=T, alternative="two.sided")$p.value,
              # KS_Pvalue = ks.test(AvgGuide[Breast], AvgGuide[!Breast],
              #                 alternative="greater")$p.value,
              BreastCN = median(CN[Breast], na.rm=T),
              OtherCN = median(CN[!Breast], na.rm=T)) %>%
    # BreastRPKM = mean(RPKM[Breast], na.rm=T),
    # OtherRPKM = mean(RPKM[!Breast], na.rm=T)) %>%
    left_join(gene_arms %>% distinct(Gene, Arm)) %>%
    mutate(Qvalue = p.adjust(Pvalue, method='fdr'))

breast_v_rest_ceres <-
    avana_gene_level %>%
    mutate(Breast = str_detect(CellLine, "_BREAST")) %>%
    group_by(Gene) %>%
    summarise(BreastMean = mean(Gene_Effect[Breast], na.rm=T),
              OtherMean = mean(Gene_Effect[!Breast], na.rm=T),
              Diff = BreastMean - OtherMean,
              Pvalue = t.test(Gene_Effect[Breast], Gene_Effect[!Breast],
                              var.equal=T, alternative="two.sided")$p.value,
              # KS_Pvalue = ks.test(Gene_Effect[Breast], Gene_Effect[!Breast],
              #                 alternative="greater")$p.value,
              BreastCN = median(CN[Breast], na.rm=T),
              OtherCN = median(CN[!Breast], na.rm=T)) %>%
    # BreastRPKM = mean(RPKM[Breast], na.rm=T),
    # OtherRPKM = mean(RPKM[!Breast], na.rm=T)) %>%
    left_join(gene_arms %>% distinct(Gene, Arm)) %>%
    mutate(Qvalue = p.adjust(Pvalue, method='fdr'))



xlimits <- c(min(breast_v_rest_avgguide$Diff, breast_v_rest_ceres$Diff, na.rm=T),
             max(breast_v_rest_avgguide$Diff, breast_v_rest_ceres$Diff, na.rm=T))
             # 0.2)


# avgguide_gene_labels <- breast_v_rest_avgguide %>%
#     filter(Qvalue < 0.01) %$% Gene
# ceres_gene_labels <- breast_v_rest_ceres %>%
#     filter(Qvalue < 0.01) %$% Gene

avgguide_gene_labels <- breast_v_rest_avgguide %>%
    filter(Qvalue < 0.01 & Diff < 0 |
               Qvalue < 0.05 & Diff < -0.3) %$% Gene
ceres_gene_labels <- breast_v_rest_ceres %>%
    filter(Qvalue < 0.01 & Diff < 0 |
               Qvalue < 0.05 & Diff < -0.25) %$% Gene


ymax <- min(breast_v_rest_avgguide$Qvalue, breast_v_rest_ceres$Qvalue)

breast_diff_dep_plot(breast_v_rest_avgguide, which_arms="chr8q", gene_labels=avgguide_gene_labels,
                     title="Uncorrected Data", xlims=xlimits,
                     break_fdr=1e-4, max_fdr=ymax) +
    ggsave(file.path(out_dir, "breast_chr8_avgguide.pdf"), width=4, height=6)
breast_diff_dep_plot(breast_v_rest_ceres, which_arms="chr8q", gene_labels=ceres_gene_labels,
                     title="CERES", xlims=xlimits,
                     break_fdr=1e-4, max_fdr=ymax) +
    ggsave(file.path(out_dir, "breast_chr8_ceres.pdf"), width=4, height=6)

# breast_diff_means_plot(breast_v_rest_avgguide, which_arms="chr8q", gene_labels=avgguide_gene_labels,
#                        title="Uncorrected Data", xlims=xlimits) +
#     ggsave(file.path(out_dir, "breast_chr8_avgguide.pdf"), width=4, height=6)
# breast_diff_means_plot(breast_v_rest_ceres, which_arms="chr8q", gene_labels=ceres_gene_labels,
#                        title="CERES", xlims=xlimits) +
#     ggsave(file.path(out_dir, "breast_chr8_ceres.pdf"), width=4, height=6)

avg_guide_sig_genes <- breast_v_rest_avgguide %>%
    filter(#Arm == "chr8q",
           Diff < 0,
           Qvalue < 0.01)

ceres_sig_genes <- breast_v_rest_ceres %>%
    filter(#Arm == "chr8q",
           Diff < 0,
           Qvalue < 0.01)

# false_positives <- setdiff(avg_guide_sig_genes$Gene, ceres_sig_genes$Gene)
# true_positives <- ceres_sig_genes$Gene


rnai_genes <- c(avg_guide_sig_genes$Gene,
                ceres_sig_genes$Gene) %>%
    c(., unlist(str_split(., "-"))) %>% unique()


rnai_validation <- rnai_dep %>%
    filter(Gene %in% rnai_genes) %>%
    mutate(Breast = str_detect(CellLine, "_BREAST")) %>%
    group_by(Gene) %>%
    summarise(BreastMean = mean(RNAi[Breast], na.rm=T),
              OtherMean = mean(RNAi[!Breast], na.rm=T),
              Diff = BreastMean - OtherMean,
              Pvalue = t.test(RNAi[Breast], RNAi[!Breast],
                              var.equal=F, alternative="two.sided")$p.value,
              # KS_Pvalue = ks.test(Gene_Effect[Breast], Gene_Effect[!Breast],
              #                 alternative="greater")$p.value,
              BreastCN = median(CN[Breast], na.rm=T),
              OtherCN = median(CN[!Breast], na.rm=T)) %>%
    left_join(gene_arms %>% distinct(Gene, GeneChr, GeneStart, GeneEnd, Arm)) %>%
    mutate(Qvalue = p.adjust(Pvalue, method='fdr'))



rnai_validation %>%
    filter(Gene %in% avg_guide_sig_genes$Gene,
           Arm == "chr8q") %>%
    mutate(Qvalue = p.adjust(Pvalue, method='fdr')) %>%
    ggplot(aes(x=Diff, y=-log10(Pvalue), color=Arm)) +
    geom_hline(yintercept=-log10(0.05), lty=2, color="grey70") +
    geom_vline(xintercept=0, lty=2, color="grey70") +
    geom_point(size=1, show.legend=F, color="#E41A1C") +
    geom_text_repel(aes(label=Gene), size=4, color="#E41A1C",
                    show.legend=F,
                    segment.color="grey70",
                    min.segment.length=unit(1, "lines"),
                    box.padding=unit(0.15, "lines")) +
    scale_y_continuous(limits=c(0, max(-log10(rnai_validation$Pvalue))*1.02),
                       expand=c(0,0)) +
    labs(x="Breast Differential Dependency (Z-score)", y="-log10(p-value)",
         title="RNAi") +
    ggsave(file.path(out_dir, "breast_chr8_rnai.pdf"), width=5, height=5)

rnai_validation %>%
    filter(Gene %in% ceres_sig_genes$Gene,
           Arm != "chr8q") %>%
    mutate(Qvalue = p.adjust(Pvalue, method='fdr')) %>%
    ggplot(aes(x=Diff, y=-log10(Pvalue), color=Arm)) +
    geom_hline(yintercept=-log10(0.05), lty=2, color="grey70") +
    geom_vline(xintercept=0, lty=2, color="grey70") +
    geom_point(size=1, show.legend=F, color="black") +
    geom_text_repel(aes(label=Gene), size=4, color="black",
                    show.legend=F,
                    segment.color="grey70",
                    min.segment.length=unit(1, "lines"),
                    box.padding=unit(0.15, "lines")) +
    scale_y_continuous(limits=c(0, max(-log10(rnai_validation$Pvalue))*1.02),
                       expand=c(0,0)) +
    labs(x="Breast Differential Dependency (Z-score)", y="-log10(p-value)",
         title="RNAi") +
    ggsave(file.path(out_dir, "breast_other_rnai.pdf"), width=5, height=5)

breast_genes <- c("TRPS1", "ESR1", "FOXA1", "GATA3",
                  "GRHL2", "TFAP2C", "SPDEF", "ZNF652")

breast_df <- avana_gene_level %>%
    filter(Gene %in% breast_genes) %>%
    group_by(Gene) %>%
    mutate(Gene_Effect = Gene_Effect - mean(Gene_Effect, na.rm=T)) %>%
    ungroup() %>%
    mutate(Gene = factor(Gene, levels=breast_genes, ordered=T),
           Breast = ifelse(str_detect(CellLine, "BREAST"), "Breast Cancer Cell Lines  ", "Other Cell Lines  ")) %>%
    arrange(desc(Breast))


ggplot(breast_df, aes(x=RPKM, y=Gene_Effect, color=Breast)) +
    facet_wrap(~ Gene, ncol=4, scales="free_x", ) +
    geom_vline(xintercept=-Inf, size=1, color="black") +
    geom_hline(yintercept=0, linetype=2, color="grey30") +
    geom_point(size=1) +
    geom_text(aes(x=X, label=Gene),
              y=max(breast_df$Gene_Effect, na.rm=T),
              color=c(brewer_pal(palette="Set1")(1), rep("black", 3)) %>%
                  rep(2),
              hjust=0, vjust=1,
              data=breast_df %>% group_by(Gene) %>% summarise(X=min(RPKM,na.rm=T))) +
    scale_color_manual(
        values=c("Other Cell Lines  "="black",
                 "Breast Cancer Cell Lines  "= brewer_pal(palette="Set1")(8)[8]),
        guide=guide_legend(override.aes=list(size=2))) +
    scale_x_continuous(breaks=seq(-2, 12, by=2)) +
    labs(x="Gene Expression (log2RPKM)", y="Differential Dependency") +
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          legend.title=element_blank(),
          legend.position="top",
          legend.margin=margin(0,0,0,0, "pt"),
          axis.line.y=element_blank(),
          axis.text=element_text(size=10)) +
    ggsave(file.path(out_dir, str_c("breast_tfs.pdf")),
           width=8, height=5)
