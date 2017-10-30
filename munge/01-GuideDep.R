
cache_helper("gecko_guide_dep", function() {
    dep <- read_csv("./data/raw/gecko_sgRNA_logFCnorm.csv") %>%
        as.data.frame() %>%
        set_rownames(.[["X1"]]) %>%
        select(-X1) %>%
        as.matrix
})

cache_helper("gecko_replicate_map", function() {
    read_tsv("./data/raw/gecko_replicate_map.tsv")
})


cache_helper("avana_guide_dep", function() {
    dep <- read_csv("./data/raw/avana_sgRNA_logFCnorm.csv") %>%
        as.data.frame() %>%
        set_rownames(.[["X1"]]) %>%
        select(-X1) %>%
        as.matrix
})

cache_helper("avana_replicate_map", function() {
    read_tsv("./data/raw/avana_replicate_map.tsv")
})
