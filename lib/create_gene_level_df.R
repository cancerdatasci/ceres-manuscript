create_gene_level_df <- function(ceres_input, ceres_output, cn_df, xpr_df,
                                 ne_genes, cce_genes,
                                 do_parallel=F) {


    guide_gene <- ceres_input$guide_locus %>%
        left_join(ceres_input$locus_gene) %>%
        filter(!is.na(Gene), Gene != "") %>%
        distinct(Guide, Gene)


    avg_guide <- ceres_input$sgRNA_effect_observed %>%
        select(Guide, CellLine, sgRNA_Observed) %>%
        df.to.mat() %>%
        collapse_rows_of_matrix(guide_gene, group_var=Gene, sample_var=Guide,
                                collapse_fun=Matrix::colMeans, na.rm=T,
                                do_parallel=do_parallel)

    gene_level <-
        make.a.tidy.dataset(
            list(Gene_Effect=df.to.mat(ceres_output$gene_effect %>%
                                           select(Gene, CellLine, Gene_Effect)),
                 AvgGuide = avg_guide,
                 CN = df.to.mat(cn_df %>%
                                    select(Gene, CellLine, CN)),
                 RPKM = df.to.mat(xpr_df %>%
                                      select(Gene, CellLine, RPKM) %>%
                                      distinct(Gene, CellLine, .keep_all=T))),
                 # MeanCuts = df.to.mat(gene_cuts %>% select(Gene, CellLine, MeanCuts))),
            dim.names=c("Gene", "CellLine"), use.dims="Gene_Effect") %>%
        filter(!is.na(CN), !is.na(AvgGuide), !is.na(Gene_Effect)) %>%
        group_by(CellLine) %>%
        mutate(AvgGuide = - (AvgGuide - median(AvgGuide[Gene %in% ne_genes])) /
                   (median(AvgGuide[Gene %in% cce_genes]) - median(AvgGuide[Gene %in% ne_genes])),
               Gene_Effect = - (Gene_Effect - median(Gene_Effect[Gene %in% ne_genes])) /
                   (median(Gene_Effect[Gene %in% cce_genes]) - median(Gene_Effect[Gene %in% ne_genes]))) %>%
        ungroup()



}
