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

load("./cache/ccds.RData")
load("./cache/cn_seg.RData")
load("./cache/cn_genes.RData")
load("./cache/lineage_df.RData")


avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")

ccds_genes <- ccds %>%
    select(Gene = gene,
           Chr = chromosome,
           Start = gene_start,
           End = gene_end) %>%
    distinct %>%
    group_by(Gene, Chr) %>%
    summarise(Start = min(Start),
              End = max(End))


min_probes <- 10
amp_thresh <- 2.8
del_thresh <- 1.2

cytobands <- read_tsv("./data/raw/hg19_cytoBand.txt",
                      col_names=c("Chr", "Start", "End", "Band", "Stain"))

acrocentric <- str_c(c(13, 14, 15, 21, 22), 'p')


arm_coords <- cytobands %>%
    mutate(ArmLabel = str_sub(Band, 1, 1)) %>%
    group_by(Chr, ArmLabel) %>%
    summarise(ArmStart = min(Start) + 1,
              ArmEnd = max(End)) %>%
    ungroup %>%
    mutate(Arm = str_c(Chr, ArmLabel))


### Make arm level calls

cn_seg_arm <-
    cn_seg %>% ldply(.id="CellLine") %>%
    left_join(lineage_df) %>%
    filter(!is.na(Lineage)) %>%
    filter(Num_Probes >= min_probes) %>%
    left_join(arm_coords, by="Chr") %>%
    mutate(End = ifelse(Start < ArmEnd & End > ArmEnd, ArmEnd, End),
           Start = ifelse(End > ArmStart & Start < ArmStart, ArmStart, Start)) %>%
    filter(Start <= ArmEnd & End >= ArmStart) %>%
    mutate(Width = (End - Start) / 1000)



arm_level_calls <- cn_seg_arm %>%
    filter(!Arm %in% acrocentric) %>%
    group_by(Lineage, CellLine, Arm) %>%
    summarise(WeightedMedian = weightedMedian(CN, Width)) %>%
    mutate(ArmCall = case_when(WeightedMedian > amp_thresh ~ 1,
                            WeightedMedian < del_thresh ~ -1,
                            TRUE ~ 0)) %>%
    ungroup


min_cell_lines <- 5

all_lineages <- arm_level_calls %>%
    distinct(Lineage, CellLine) %>%
    group_by(Lineage) %>%
    summarise(n = sum(unique(CellLine) %in% unique(avana_gene_level$CellLine))) %>%
    filter(n >= min_cell_lines) %$%
    Lineage %>% set_names(., .)


lineage_arm_cn_results <- ldply(all_lineages, function(cur_lin) {
    # print(cur_lin)
    lin_df <- arm_level_calls %>%
        mutate(is_cur_lin = Lineage == cur_lin)
    res <- lin_df %>%
        group_by(Arm) %>%
        filter(any(ArmCall == 1)) %>%
        do(fisher.test(.$is_cur_lin, .$ArmCall == 1,
                       alternative = 'greater') %>%
               broom::tidy())
}, .id = 'Lineage') %>%
    mutate(q.value = p.adjust(p.value, method = 'fdr'))



lineage_v_rest <- all_lineages %>%
    ldply(function(lineage) {
        print(lineage)
        lineage_lines <- lineage_df %>% filter(Lineage == lineage) %$% CellLine

        avana_before_ceres <- avana_gene_level %>%
            mutate(Lineage = CellLine %in% lineage_lines) %>%
            group_by(Gene) %>%
            do(broom::glance(t.test(.$AvgGuide[.$Lineage], .$AvgGuide[!.$Lineage],
                                    alternative="less", var.equal=T))) %>%
            mutate(Dataset="Uncorrected Data")

        avana_after_ceres <- avana_gene_level %>%
            mutate(Lineage = CellLine %in% lineage_lines) %>%
            group_by(Gene) %>%
            do(broom::glance(t.test(.$Gene_Effect[.$Lineage], .$Gene_Effect[!.$Lineage],
                                    alternative="less", var.equal=T))) %>%
            mutate(Dataset="CERES")

        lin_result <- bind_rows(avana_before_ceres, avana_after_ceres)
    }, .id="Lineage", .parallel=do_parallel)



lineage_v_rest_arm <- lineage_v_rest %>%
    group_by(Lineage, Dataset) %>%
    mutate(q.value = p.adjust(p.value, method='fdr')) %>%
    ungroup %>%
    left_join(ccds_genes, by="Gene") %>%
    left_join(arm_coords, by="Chr") %>%
    filter(Start < ArmEnd, End > ArmStart) %>%
    select(-ArmLabel, -ArmStart, -ArmEnd, -Chr, -Start, -End)



arm_cn_sig_threshold <- 0.05
enriched_arms <- lineage_arm_cn_results %>%
    filter(q.value < arm_cn_sig_threshold) %>%
    arrange(desc(estimate))

diff_dep_sig_threshold <- 0.05
enrichment_results <-
    enriched_arms %>%
    ddply(.(Lineage, Arm), function(enriched_arm) {
        test_lineage <- as.character(enriched_arm$Lineage[1])
        test_arm <- as.character(enriched_arm$Arm[1])
        lineage_v_rest_arm %>%
            filter(Lineage == test_lineage) %>%
            mutate(SigDep = q.value < diff_dep_sig_threshold &
                       estimate1 < estimate2,
                   OnArm = Arm == test_arm) %>%
            group_by(Dataset) %>%
            summarise(ArmDeps = sum(SigDep & OnArm),
                      OtherDeps = sum(SigDep & !OnArm),
                      ArmNondep = sum(!SigDep & OnArm),
                      OtherNondep = sum(!SigDep & !OnArm))
    }) %>%
    mutate(TotalDeps = ArmDeps + OtherDeps,
           FractionDeps = ArmDeps/TotalDeps,
           LinArm = paste(Lineage, Arm))



enrichment_results %>%
    group_by(Lineage, Arm) %>%
    filter(any(TotalDeps >= 10)) %>%
    mutate(LinArm = factor(LinArm,
                           levels=filter(., Dataset=="Uncorrected Data") %>%
                               arrange(desc(FractionDeps)) %$% LinArm, ordered=T),
           Dataset = factor(Dataset, levels=c("Uncorrected Data", "CERES", ordered=T))) %>%
    ggplot(aes(x=LinArm, y=FractionDeps, color=Dataset)) +
    geom_point(size=2) +
    scale_color_manual(values=brewer_pal(palette="Set1")(2)) +
    scale_y_continuous(limits=c(0, 0.4), expand=c(0,0), labels=scales::percent) +
    labs(y = "Differential Dependencies\non Chromosome Arm (%)") +
    theme(panel.grid.major.x=element_line(color="grey80", linetype=1, size=0.25),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1),
          axis.text.y=element_text(size=10),
          legend.position=c(0.5, 1),
          legend.justification=c(0.5, 1),
          legend.background=element_rect(fill="white", linetype=0),
          legend.title=element_blank()) +
    ggsave(file.path(out_dir, "lineage_arm_fraction_avana.pdf"), width=8, height=4)



fisher_results <- enrichment_results %>%
    group_by(LinArm, Dataset) %>%
    do(broom::glance(fisher.test(matrix(c(.$ArmDeps, .$OtherDeps,
                                          .$ArmNondep, .$OtherNondep), nrow=2)))) %>%
    left_join(enrichment_results) %>%
    group_by(Lineage, Arm) %>%
    filter(any(TotalDeps >= 10)) %>%
    mutate(LinArm = factor(LinArm,
                           levels=filter(., Dataset=="Uncorrected Data") %>%
                               arrange(desc(estimate)) %$% LinArm, ordered=T),
           Dataset = factor(Dataset, levels=c("Uncorrected Data", "CERES", ordered=T))) %>%
    mutate(X = as.numeric(LinArm) + ifelse(Dataset=="Uncorrected Data", -0.1, 0.1))

ggplot(fisher_results, aes(x=X, y=estimate, color=Dataset)) +
    geom_segment(aes(x=X, xend=X,
                     y=conf.low, yend=conf.high)) +
    geom_point(size=2) +
    scale_color_manual(values=brewer_pal(palette="Set1")(2)) +
    scale_y_continuous(limits=c(0, max(fisher_res_plot$conf.high)), expand=c(0,0)) +
    scale_x_continuous(breaks=1:length(levels(fisher_res_plot$LinArm)),
                       labels=levels(fisher_res_plot$LinArm),
                       expand=c(0, 0.5)) +
    labs(y = "Enrichment of Dependencies\non Chromosome Arm") +
    theme(panel.grid.major.x=element_line(color="grey80", linetype=1, size=0.25),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=12),
          axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1),
          axis.text.y=element_text(size=10),
          legend.position=c(0.5, 1),
          legend.justification=c(0.5, 1),
          legend.background=element_rect(fill="white", linetype=0),
          legend.title=element_blank()) +
    ggsave(file.path(out_dir, "lineage_arm_enrichment_avana.pdf"), width=8, height=4)

