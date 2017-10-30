make_ceres_inputs <- function(inputs_dir, guide_dep_mat, neg_ctrl_guides=NULL, pos_ctrl_guides=NULL,
                              guide_cn_mat, guide_alns_df,  gene_df,
                              rep_map=NULL, collapse_replicates=F,
                              remove_XY=F,
                              cell_lines=NULL, dep_normalize="none",
                              remove_lines=NULL) {


    dir.create(inputs_dir, recursive=T, showWarnings=F)

    guide_dep <- guide_dep_mat %>%
        set_rownames(str_extract(rownames(.), "^[ACGT]+")) %>%
        {.[unique(rownames(.)),]} %>%
        remove.rows.all.nas()



    locus_cn <- guide_cn_mat[str_detect(rownames(guide_cn_mat), "chr"), ] %>%
        set_rownames(rownames(.) %>% str_extract("chr.+$")) %>%
        {.[rownames(.),]} %>%
        remove.rows.all.nas()

    non_targeting_cn <- guide_cn_mat[!str_detect(rownames(guide_cn_mat), "chr"), ]

    guide_locus_df <- guide_alns_df %>%
        transmute(Guide,
                  Locus = str_c(rname, Cut.Pos, strand, sep="_"),
                  Value = 1) %>%
        distinct()

    locus_gene_df <- gene_df %>%
        transmute(Locus = str_c(Chr, Cut.Pos, Strand, sep="_"),
                  Gene,
                  Value = 1) %>%
        distinct()

    if (is.null(rep_map)) {
        cell_lines_to_use <- intersect(colnames(guide_dep), colnames(locus_cn)) %>%
            setdiff(remove_lines)
    } else {
        cell_lines_to_use <- intersect(rep_map$CellLine, colnames(locus_cn)) %>%
            setdiff(remove_lines)
    }

    if (!is.null(cell_lines)) {
        cell_lines_to_use <- intersect(cell_lines, cell_lines_to_use)
    }


    loci_to_use <- intersect(guide_locus_df$Locus, rownames(locus_cn))
    guides_to_use <-
        intersect(rownames(guide_dep),
                  c(guide_locus_df %>% filter(Locus %in% loci_to_use) %$% Guide,
                    rownames(non_targeting_cn)))


    if (remove_XY) {
        XY_guides <- gene_df %>%
            group_by(Guide) %>%
            filter(any(Chr == "chrX" | Chr == "chrY"))

        guides_to_use <- setdiff(guides_to_use, XY_guides$Guide)
    }


    if (is.null(rep_map)) {
        # if (!is.null(neg_ctrl_guides)) {
        #     negctrl_dep <- negctrl_dep[,cell_lines_to_use]
        # }
        guide_dep <- guide_dep[guides_to_use, cell_lines_to_use]

    } else {
        # if (!is.null(neg_ctrl_guides)) {
        #
        #     negctrl_dep <- negctrl_dep[,rep_map %>%
        #                                    filter(CellLine %in% cell_lines_to_use) %$%
        #                                    Replicate]
        # }
        guide_dep <- guide_dep[guides_to_use, rep_map %>%
                                   filter(CellLine %in% cell_lines_to_use) %$%
                                   Replicate]
    }

    if (dep_normalize=="posneg") {
        guide_dep <- guide_dep %>%
            as.data.frame %>%
            mutate(Guide = rownames(.)) %>%
            gather(Sample, Dep, -Guide) %>%
            group_by(Sample) %>%
            mutate(NegCtrlMedian = median(Dep[Guide %in% neg_ctrl_guides], na.rm=T),
                   PosCtrlMedian = median(Dep[Guide %in% pos_ctrl_guides], na.rm=T)) %>%
            mutate(Dep = - (Dep - NegCtrlMedian) / (PosCtrlMedian - NegCtrlMedian)) %>%
            ungroup %>%
            select(Sample, Guide, Dep) %>%
            spread(Sample, Dep) %>%
            as.data.frame %>%
            set_rownames(.$Guide) %>%
            select(-Guide) %>%
            as.matrix
    } else if (dep_normalize=="zmad") {
        guide_dep <- guide_dep %>%
            as.data.frame %>%
            mutate(Guide = rownames(.)) %>%
            gather(Sample, Dep, -Guide) %>%
            group_by(Sample) %>%
            mutate(Dep = (Dep - median(Dep, na.rm=T)) / mad(Dep, na.rm=T)) %>%
            ungroup %>%
            select(Sample, Guide, Dep) %>%
            spread(Sample, Dep) %>%
            as.data.frame %>%
            set_rownames(.$Guide) %>%
            select(-Guide) %>%
            as.matrix
    } else if (dep_normalize=="none") {

    } else {
        stop("Error: normalization not recognized")
    }


    if (collapse_replicates) {
        if (is.null(rep_map)) stop("Error: no replicate map")
        guide_dep <- collapse_reps(guide_dep, rep_map)
    }


    guide_locus_df <- guide_locus_df %>%
        filter(Guide %in% guides_to_use,
               Locus %in% loci_to_use)
    locus_gene_df <- locus_gene_df %>%
        filter(Locus %in% loci_to_use)

    gene_ids <- gene_df %>%
        select(GeneID,
               Gene) %>%
        distinct()

    saveRDS(guide_dep, file.path(inputs_dir, "guide_sample_dep.Rds"))
    # if (!is.null(neg_ctrl_guides)) {
    #     saveRDS(negctrl_dep, file.path(inputs_dir, "negctrl_guide_sample_dep.Rds"))
    # }
    saveRDS(guide_locus_df, file.path(inputs_dir, "guide_locus_map3.Rds"))
    saveRDS(locus_gene_df, file.path(inputs_dir, "locus_gene_map3.Rds"))
    saveRDS(locus_cn, file.path(inputs_dir, "locus_sample_cn.Rds"))
    saveRDS(gene_ids, file.path(inputs_dir, "gene_id_table.Rds"))

    if (!is.null(rep_map) & !collapse_replicates) {
        rep_map %>% filter(Replicate %in% colnames(guide_dep)) %>%
            saveRDS(file.path(inputs_dir, "replicate_map.Rds"))
        collapse_reps(guide_dep, rep_map) %>%
            saveRDS(file.path(inputs_dir, "guide_dep_reps_collapsed.Rds"))
    }


    invisible(NULL)
}




get_gene_annotations <- function(genes, guide_alns,
                                 chromosomes=paste0("chr",c(as.character(1:22), "X", "Y"))) {


    gene_annot_grs <- genes %>%
        makeGRangesFromDataFrame(seqinfo=hg19info, keep.extra.columns=T)



    guide_aln_grs <- guide_alns %>%
        dplyr::select(Guide, Chr=rname, Start=Cut.Pos, strand) %>%
        mutate(End = Start) %>%
        filter(Chr %in% chromosomes) %>%
        distinct() %>%
        makeGRangesFromDataFrame(seqinfo=hg19info, keep.extra.columns=T)

    hits <- findOverlaps(guide_aln_grs, gene_annot_grs, ignore.strand=T) %>%
        as.data.frame

    gene_df <- hits %>%
        transmute(Guide = guide_aln_grs$Guide[queryHits],
                  Chr = seqnames(guide_aln_grs)[queryHits] %>% as.character(),
                  Cut.Pos = start(guide_aln_grs)[queryHits] %>% as.integer(),
                  Strand = strand(guide_aln_grs)[queryHits] %>% as.character(),
                  Gene = gene_annot_grs$gene[subjectHits],
                  GeneID = gene_annot_grs$gene_id[subjectHits],
                  CDS_Strand = strand(gene_annot_grs)[subjectHits] %>% as.character(),
                  CDS_Start = start(gene_annot_grs)[subjectHits] %>% as.integer(),
                  CDS_End = end(gene_annot_grs)[subjectHits] %>% as.integer()) %>%
        distinct
}


collapse_reps <- function(replicate_dep, replicate_map, do_parallel=F) {
    rep_map <- replicate_map %>% filter(Replicate %in% colnames(replicate_dep))
    collapsed <- daply(rep_map, .(CellLine),
                       function(reps) {
                           rowMeans(replicate_dep[,reps$Replicate, drop=F], na.rm=T)
                       }, .parallel=do_parallel) %>% t()
}

