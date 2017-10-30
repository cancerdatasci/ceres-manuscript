threads <- 2

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- T
} else {
    do_parallel <- F
}


cache_helper("gecko_guide_cn", function() {

    cn_seg <- cn_seg[names(cn_seg) %in% gecko_replicate_map$CellLine]

    cn_seg_gr <- llply(cn_seg, makeGRangesFromDataFrame,
                       seqinfo=hg19info, keep.extra.columns=T)

    gecko_guide_alns_gr <- gecko_guide_alns %>%
        filter(!is.na(rname)) %>%
        mutate(Chr = rname,
               Start = Cut.Pos,
               End = Cut.Pos,
               AlnID = str_c(Guide, Chr, Start, strand, sep="_")) %>%
        distinct(Guide, Chr, Start, End, AlnID) %>%
        makeGRangesFromDataFrame(seqinfo=hg19info, keep.extra.columns=T)

    gecko_guide_no_alns <- gecko_guide_alns %>%
        filter(is.na(rname)) %>%
        distinct(Guide, .keep_all=T)

    gecko_guide_cn <-
        intersect_guide_with_copy_number(gecko_guide_alns_gr, cn_seg_gr,
                                         CN.column="CN",
                                         guide.column="AlnID",
                                         do_parallel=do_parallel) %>%
        rbind(matrix(0, dimnames=list(gecko_guide_no_alns$Guide, colnames(.)),
                     nrow=nrow(gecko_guide_no_alns), ncol=ncol(.)))

})

cache_helper("avana_guide_cn", function() {

    cn_seg <- cn_seg[names(cn_seg) %in% avana_replicate_map$CellLine]

    cn_seg_gr <- llply(cn_seg, makeGRangesFromDataFrame,
                       seqinfo=hg19info, keep.extra.columns=T)

    avana_guide_alns_gr <- avana_guide_alns %>%
        filter(!is.na(rname)) %>%
        mutate(Chr = rname,
               Start = Cut.Pos,
               End = Cut.Pos,
               AlnID = str_c(Guide, Chr, Start, strand, sep="_")) %>%
        distinct(Guide, Chr, Start, End, AlnID) %>%
        makeGRangesFromDataFrame(seqinfo=hg19info, keep.extra.columns=T)

    avana_guide_no_alns <- avana_guide_alns %>%
        filter(is.na(rname)) %>%
        distinct(Guide, .keep_all=T)

    avana_guide_cn <-
        intersect_guide_with_copy_number(avana_guide_alns_gr, cn_seg_gr,
                                         CN.column="CN",
                                         guide.column="AlnID",
                                         do_parallel=do_parallel) %>%
        rbind(matrix(0, dimnames=list(avana_guide_no_alns$Guide, colnames(.)),
                     nrow=nrow(avana_guide_no_alns), ncol=ncol(.)))

})



registerDoSEQ()
