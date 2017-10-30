library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))

threads <- 2

if (threads > 1) {
    library(doMC)
    registerDoMC(cores=threads)
    do_parallel <- TRUE
} else {
    do_parallel <- FALSE
}

out_dir <- file.path("./output/figures/example_regions", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


load("./cache/hg19info.RData")
load("./cache/cn_seg.RData")
load("./cache/ccds.RData")

avana_gene_level <- readRDS("./data/ceres_cache_final/avana_gene_level.rds")
gecko_gene_level <- readRDS("./data/ceres_cache_final/gecko_gene_level.rds")
wang_gene_level <- readRDS("./data/ceres_cache_final/wang_gene_level.rds")

ccds_genes <- ccds %>%
    select(Gene = gene,
           Chr = chromosome,
           # Strand = strand,
           Start = gene_start,
           End = gene_end) %>%
    distinct %>%
    group_by(Gene, Chr) %>%
    summarise(Start = min(Start),
              End = max(End))


min_probes <- 10
cn_lim <- c(0, 15)
crispr_lim <- c(-2, 1)


plot_regions <- read_tsv("./data/manual/example_regions.tsv") %>%
    mutate(Selector = str_c(AmpName, Sample, sep=" "))


cn <- ldply(cn_seg, .id="CellLine") %>%
    filter(Num_Probes > min_probes)


cn_gr_list <- dlply(cn, .(CellLine), makeGRangesFromDataFrame,
                    keep.extra.columns=T, seqinfo=hg19info)



for (i in 1:nrow(plot_regions)) {
    plot_region <- plot_regions[i,] %>%
        makeGRangesFromDataFrame(seqinfo=hg19info, keep.extra.columns=T)

    cell_line <- plot_region$Sample

    # if (!cell_line %in% unique(avana_gene_level$CellLine)) next

    if (cell_line %in% unique(avana_gene_level$CellLine)) {
        gene_level <- avana_gene_level %>% filter(CellLine == cell_line)
        plot_dir <- file.path(out_dir, "example_regions_avana")
    } else if (cell_line %in% unique(gecko_gene_level$CellLine)) {
        gene_level <- gecko_gene_level %>% filter(CellLine == cell_line)
        plot_dir <- file.path(out_dir, "example_regions_gecko")
    } else if (cell_line %in% unique(wang_gene_level$CellLine)) {
        gene_level <- wang_gene_level %>% filter(CellLine == cell_line)
        plot_dir <- file.path(out_dir, "example_regions_wang")
    } else {
        next
    }

    before_ceres <- gene_level %>%
        mutate(Dep = AvgGuide) %>%
        filter(CellLine == cell_line,
               !is.na(Gene), Gene != "") %>%
        left_join(ccds_genes, by="Gene")

    after_ceres <- gene_level %>%
        mutate(Dep = Gene_Effect) %>%
        filter(CellLine == cell_line,
               !is.na(Gene), Gene != "") %>%
        left_join(ccds_genes, by="Gene")

    before_gr <- makeGRangesFromDataFrame(before_ceres, keep.extra.columns=T, seqinfo=hg19info)
    after_gr <- makeGRangesFromDataFrame(after_ceres, keep.extra.columns=T, seqinfo=hg19info)
    cn_gr <- cn_gr_list[[cell_line]]

    cn_gr$Before_Median <- intersect_data_with_segment(cn_gr, before_gr, dat.column="Dep",
                                                    func=median, na.value=NA, na.rm=T)
    cn_gr$Before_Num <- intersect_data_with_segment(cn_gr, before_gr, dat.column="Dep",
                                                 func=length, na.value=NA)
    cn_gr$After_Median <- intersect_data_with_segment(cn_gr, after_gr, dat.column="Dep",
                                                   func=median, na.value=NA, na.rm=T)
    cn_gr$After_Num <- intersect_data_with_segment(cn_gr, after_gr, dat.column="Dep",
                                                func=length, na.value=NA)


    before_gviz <- cn_dep_gviz(plot_region,
                               crispr=before_gr, rnai=NULL, cn=cn_gr,
                               highlight.gene=plot_region$Gene,
                               crispr.seg.column="Before_Median",
                               sample.name=plot_region$CellLine,
                               chr.name=plot_region$AmpName,
                               CN.column="CN",
                               CN.baseline=2,
                               CN.track.title="Copy Number",
                               CRISPR.track.title="Uncorrected",
                               close.gaps=TRUE,
                               close.NA.segs=TRUE,
                               CRISPR.column="Dep",
                               CN.ylim=cn_lim,
                               CRISPR.ylim=crispr_lim,
                               CRISPR.seg.lwd=2,
                               highlight.gene.cex=1.25,
                               cex.main=1.25,
                               cex.title=1.25, cex.axis=1,
                               innerMargin=20,
                               fontface.title=1,
                               fontfamily="Helvetica",
                               add_gridlines=F)


    after_gviz <- cn_dep_gviz(plot_region,
                              crispr=after_gr, rnai=NULL, cn=cn_gr,
                              highlight.gene=plot_region$Gene,
                              crispr.seg.column="After_Median",
                              sample.name=plot_region$CellLine,
                              chr.name=plot_region$AmpName,
                              CN.column="CN",
                              CN.baseline=2,
                              CN.track.title="Copy Number",
                              CRISPR.track.title="CERES",
                              close.gaps=TRUE,
                              close.NA.segs=TRUE,
                              CRISPR.column="Dep",
                              CN.ylim=cn_lim,
                              CRISPR.ylim=crispr_lim,
                              CRISPR.seg.lwd=2,
                              highlight.gene.cex=1.25,
                              cex.main=1.25,
                              cex.title=1.25, cex.axis=1,
                              innerMargin=20,
                              fontface.title=1,
                              fontfamily="Helvetica",
                              add_gridlines=F)


    dir.create(plot_dir, recursive=T, showWarnings=F)

    png(file.path(plot_dir,
                  str_c(str_replace(plot_region$Selector, "\\s+", "-"), ".png")),
        width=960, height=480)

    pushViewport(viewport(layout=grid.layout(1, 2)))
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=1))

    do.call(Gviz::plotTracks, c(before_gviz, add=T))
    popViewport()
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))

    do.call(Gviz::plotTracks, c(after_gviz, add=T))
    popViewport()

    dev.off()

    pdf(file.path(plot_dir,
                  str_c(str_replace(plot_region$Selector, "\\s+", "-"), "-avg_guide.pdf")),
        width=5, height=4)
    do.call(Gviz::plotTracks, c(before_gviz, add=T))
    dev.off()

    pdf(file.path(plot_dir,
                  str_c(str_replace(plot_region$Selector, "\\s+", "-"), "-ceres.pdf")),
        width=5, height=4)
    do.call(Gviz::plotTracks, c(after_gviz, add=T))
    dev.off()


}

