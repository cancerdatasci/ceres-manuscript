
# Map guides to genome

cache_helper("gecko_guide_alns", function() {

    rownames(gecko_guide_dep) %>% unique %>% set_names(.,.) %>%
        Biostrings::DNAStringSet() %>%
        Biostrings::writeXStringSet("data/interim/gecko_guides.fa")

    system("bowtie -t -p 4 -a -v 0 -f -S hg19 data/interim/gecko_guides.fa data/interim/gecko-v0.sam")
    system("samtools view -bS -o data/interim/gecko-v0.bam data/interim/gecko-v0.sam")

    alns <- guideAlignments("data/interim/gecko-v0.bam", max.alns=100,
                    include.no.align=T, as.df=T)
})


cache_helper("avana_guide_alns", function() {
    rownames(avana_guide_dep) %>% unique %>% set_names(.,.) %>%
        Biostrings::DNAStringSet() %>%
        Biostrings::writeXStringSet("data/interim/avana_guides.fa")

    system("bowtie -t -p 4 -a -v 0 -f -S hg19 data/interim/avana_guides.fa data/interim/avana-v0.sam")
    system("samtools view -bS -o data/interim/avana-v0.bam data/interim/avana-v0.sam")

    guide_bam <- "data/interim/avana-v0.bam"
    alns <- guideAlignments("data/interim/avana-v0.bam", max.alns=100,
                            include.no.align=T, as.df=T)
})


cache_helper("gecko_genes", function() {
    get_gene_annotations(ccds, gecko_guide_alns)
})

cache_helper("avana_genes", function() {
    get_gene_annotations(ccds, avana_guide_alns)
})
