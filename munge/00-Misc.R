chromosomes <- c(as.character(1:22), "X", "Y") %>% str_c("chr", .)

cache_helper("hg19info", function() {
    Seqinfo(genome="hg19")[chromosomes]
})

cache_helper("ccds", function() {
    ccds_h37 <- read_tsv("./data/raw/CCDS.current.txt",
                         col_types=cols("#chromosome" = col_character(),
                                        "cds_from" = col_integer(),
                                        "cds_to" = col_integer())) %>%
        dplyr::rename(chromosome=`#chromosome`) %>%
        mutate(chromosome = str_c("chr", chromosome)) %>%
        filter(ccds_status %in% c("Public", "Reviewed, update pending", "Under review, update"),
               chromosome %in% chromosomes,
               !is.na(cds_from), !is.na(cds_to))

    ccds_exon_h37 <-  ccds_h37 %>%
        mutate(cds_interval = str_replace_all(cds_locations, "[\\[\\]]", "") %>%
                   str_split("\\s*,\\s*")) %>%
        unnest(cds_interval) %>%
        group_by(gene, gene_id, cds_locations) %>%
        mutate(exon_code = ifelse(cds_strand=="+", 1:n(), n():1)) %>%
        ungroup %>%
        mutate(cds_start = str_extract(cds_interval, "^[0-9]+") %>% as.integer,
               cds_end = str_extract(cds_interval, "[0-9]+$") %>% as.integer) %>%
        select(gene, gene_id, chromosome, start=cds_start, end=cds_end, strand=cds_strand,
               gene_start=cds_from, gene_end=cds_to, exon_code)
})




cache_helper("cn_seg", function() {
    read_tsv("./data/raw/CCLE_copynumber_2013-12-03.seg.txt",
             col_types="ccddid") %>%
        dplyr::rename(Chr = Chromosome) %>%
        mutate(Chr = str_c("chr", Chr),
               Start = as.integer(Start),
               End = as.integer(End)) %>%
        mutate(CN = 2 * 2^Segment_Mean) %>%
        filter(!str_detect(Chr, "chr[XY]")) %>%
        dlply(.(CCLE_name))
})



# cache_helper("ccle_mut", function() {
#     read_tsv("./data/raw/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf") %>%
#         select(CellLine = Tumor_Sample_Barcode,
#                Gene = Hugo_Symbol,
#                Chromosome, Start_position, End_position, Strand,
#                Variant_Classification, Variant_Type,
#                Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
#                Reference_allele_reads_count, Alternative_allele_reads_count)
# })



cache_helper("ccle_mut", function() {
    # ccle_mut <- load.from.taiga(data.name="ccle-mutation-data",
    #                             data.version=1) %>%
    read_tsv("./data/raw/ccle2maf_081117.txt",
             col_types = cols(Tumor_Sample_Barcode=col_character(),
                              Hugo_Symbol = col_character(),
                              Entrez_Gene_Id = col_character(),
                              Variant_Classification = col_character(),
                              Variant_Type = col_character(),
                              cDNA_Change = col_character(),
                              Protein_Change = col_character(),
                              isDeleterious = col_logical(),
                              .default = col_skip())) %>%
    dplyr::rename(CellLine = Tumor_Sample_Barcode,
                  Gene = Hugo_Symbol)
})

cache_helper("xpr", function() {
    read.gct("./data/raw/CCLE_RNAseq_081117.rpkm.gct") %>%
        set_rownames(attr(., "Description")) %>%
        {.[rownames(.) %in% ccds$gene,]} %>%
        log2 %>%
        mat.to.df("Gene", "CellLine", "RPKM")
})


cache_helper("cce_genes", function() {
    hart_ctrls <- read_tsv("./data/manual/hart2014_seed_core_essentials.txt")
    cce_genes <- hart_ctrls$`ConstitutiveCoreEssentials(CCE)` %>% str_subset("[A-Z0-9]+")
})

cache_helper("ne_genes", function() {
    hart_ctrls <- read_tsv("./data/manual/hart2014_seed_core_essentials.txt")
    ne_genes <- hart_ctrls$`Nonessential Genes (NE)` %>% str_subset("[A-Z0-9]+")
})

cache_helper("lineage_df", function() {
    lineage_pretty <- read_tsv("./data/manual/lineage_pretty.tsv")
    read_tsv("./data/manual/ccle_lineages.tsv") %>%
        dplyr::rename(Lineage = lineage, CellLine = CCLE_ID) %>%
        left_join(lineage_pretty) %>%
        mutate(Lineage = ifelse(is.na(LineagePretty), Lineage, LineagePretty))
})

