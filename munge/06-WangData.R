
threads <- 2

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- T
} else {
    do_parallel <- F
}


cache_helper("wang_guide_dep", function() {


    cell_line_labels <- read_tsv("./data/raw/wang_cell_line_labels.txt")

    sgrna_annotations <- read_tsv("./data/raw/wang_sgRNA_annotations.txt", skip=1) %>%
        select(GuideID = `sgRNA ID`,
               Guide = `sgRNA sequence`)

    normalized_counts <- read_tsv("./data/raw/wang_pool.normalized.counts.txt",
                                  col_types=cols(.default = col_number(),
                                                 sgRNA = col_character()))

    normalized_counts_tidy <- normalized_counts %>%
        gather(Expt, Count, -sgRNA)


    initial_reference <- normalized_counts_tidy %>%
        filter(Expt %in% cell_line_labels$Initial) %>%
        group_by(sgRNA) %>%
        summarise(InitialCount = sum(Count)) %>%
        ungroup %>%
        mutate(InitialCount = (InitialCount+1) * 1e6 / sum(InitialCount))


    log_fc <- normalized_counts_tidy %>%
        filter(Expt %in% cell_line_labels$Final) %>%
        left_join(cell_line_labels, by=c("Expt"="Final")) %>%
        group_by(Replicate) %>%
        mutate(FinalCount = (Count + 1) * 1e6 / sum(Count)) %>%
        ungroup %>%
        left_join(initial_reference) %>%
        mutate(logFC = log2((FinalCount)/(InitialCount)))

    guide_dep <- log_fc %>%
        dplyr::rename(GuideID=sgRNA) %>%
        left_join(sgrna_annotations) %>%
        group_by(CellLine, Guide) %>%
        summarise(logFC = mean(logFC, na.rm=T)) %>%
        ungroup() %>% select(Guide, CellLine, logFC) %>% df.to.mat()
})


cache_helper("wang_guide_alns", function() {
    rownames(wang_guide_dep) %>% unique %>% set_names(.,.) %>%
        Biostrings::DNAStringSet() %>%
        Biostrings::writeXStringSet("data/interim/wang_guides.fa")

    system("bowtie -t -p 4 -a -v 0 -f -S hg19 data/interim/wang_guides.fa data/interim/wang-v0.sam")
    system("samtools view -bS -o data/interim/wang-v0.bam data/interim/wang-v0.sam")

    wang_guide_alns <- guideAlignments("data/interim/wang-v0.bam", max.alns=100,
                                       include.no.align=T, as.df=T)
})


cache_helper("wang_genes", function() {
    wang_genes <- get_gene_annotations(ccds, wang_guide_alns)
})

cache_helper("wang_neg_ctrl_guides", function() {
    wang_neg_ctrl_guides <- wang_guide_alns %>% filter(is.na(rname)) %>%
        select(Guide, Gene)
})


cache_helper("wang_pos_ctrl_guides", function() {

    essential_gene_sets <- c("KEGG_RIBOSOME", "KEGG_PROTEASOME", "KEGG_SPLICEOSOME")

    c2 <- read.gmt("./data/raw/c2.all.v6.0.symbols.gmt", as.df=F)

    pos.ctrl.genes <- c2[essential_gene_sets] %>% unlist


    wang_pos_ctrl_guides <- wang_genes %>%
        select(Guide, Gene) %>%
        filter(Gene %in% pos.ctrl.genes) %>%
        select(Guide, Gene) %>% distinct()
})



cache_helper("wang_guide_cn", function() {
    cn_seg_gr <- cn_seg[names(cn_seg) %in% colnames(wang_guide_dep)] %>%
        llply(makeGRangesFromDataFrame,
              seqinfo=hg19info, keep.extra.columns=T)

    wang_guide_alns_gr <- wang_guide_alns %>%
        filter(!is.na(rname)) %>%
        mutate(Chr = rname,
               Start = Cut.Pos,
               End = Cut.Pos,
               AlnID = str_c(Guide, Chr, Start, strand, sep="_")) %>%
        distinct(Guide, Chr, Start, End, AlnID) %>%
        makeGRangesFromDataFrame(seqinfo=hg19info, keep.extra.columns=T)

    wang_guide_no_alns <- wang_guide_alns %>%
        filter(is.na(rname)) %>%
        distinct(Guide, .keep_all=T)

    wang_guide_cn <-
        intersect_guide_with_copy_number(wang_guide_alns_gr, cn_seg_gr,
                                         CN.column="CN",
                                         guide.column="AlnID",
                                         do_parallel=do_parallel) %>%
        rbind(matrix(0, dimnames=list(wang_guide_no_alns$Guide, colnames(.)),
                     nrow=nrow(wang_guide_no_alns), ncol=ncol(.)))

})


registerDoSEQ()

