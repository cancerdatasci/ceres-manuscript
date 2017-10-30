cache_helper("gecko_neg_ctrl_guides", function() {
    neg_ctrl_guides <- gecko_guide_alns %>% filter(is.na(rname)) %>%
        select(Guide, Gene)
})


cache_helper("avana_neg_ctrl_guides", function() {
    neg_ctrl_guides <- avana_guide_alns %>% filter(is.na(rname)) %>%
        select(Guide, Gene)
})


cache_helper("gecko_pos_ctrl_guides", function() {

    essential_gene_sets <- c("KEGG_RIBOSOME", "KEGG_PROTEASOME", "KEGG_SPLICEOSOME")

    c2 <- read.gmt("./data/raw/c2.all.v6.0.symbols.gmt", as.df=F)

    pos.ctrl.genes <- c2[essential_gene_sets] %>% unlist

    pos_ctrl_guides <- gecko_genes %>%
        select(Guide, Gene) %>%
        filter(Gene %in% pos.ctrl.genes) %>%
        select(Guide, Gene) %>% distinct()
})



cache_helper("avana_pos_ctrl_guides", function() {

    essential_gene_sets <- c("KEGG_RIBOSOME", "KEGG_PROTEASOME", "KEGG_SPLICEOSOME")

    c2 <- read.gmt("./data/raw/c2.all.v6.0.symbols.gmt", as.df=F)

    pos.ctrl.genes <- c2[essential_gene_sets] %>% unlist

    pos_ctrl_guides <- gecko_genes %>%
        select(Guide, Gene) %>%
        filter(Gene %in% pos.ctrl.genes) %>%
        select(Guide, Gene) %>% distinct()

    pos_ctrl_guides <- avana_genes %>%
        select(Guide, Gene) %>%
        filter(Gene %in% pos.ctrl.genes) %>%
        select(Guide, Gene) %>% distinct()
})
