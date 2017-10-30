threads <- 2

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- T
} else {
    do_parallel <- F
}



cache_helper("cn_genes", function() {

    gene_gr <- ccds %>%
        distinct(gene, chromosome, strand, gene_start, gene_end) %>%
        group_by(gene, chromosome, strand) %>%
        summarise(start = min(gene_start),
                  end = max(gene_end)) %>%
        makeGRangesFromDataFrame(keep.extra.columns=T)
    seg_gr <- cn_seg %>%
        llply(makeGRangesFromDataFrame, keep.extra.columns=T, .parallel=do_parallel)
    cn_mat <- intersect_gene_with_copy_number(gene_gr, seg_gr,
                                                   CN.column="CN",
                                                   gene.column="gene",
                                                   do_parallel=do_parallel) %>%
                                                   {.[!duplicated(rownames(.)),]} %>%
        remove.rows.all.nas()


    cn_genes <- mat.to.df(cn_mat, "Gene", "CellLine", "CN")

})

registerDoSEQ()
