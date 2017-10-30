
#' Plot copy number and dependency data on genomic coordinates
#'
#' @param plot.region GRanges object with one interval
#' @param crispr GRanges object with dependency data column
#' @param rnai GRanges object with dependency data column
#' @param cn GRanges object with copy number data column and optional dependency segment columns
#' @param highlight.gene gene name of single gene to highlight on plots
#' @param crispr.seg.column
#' @param rnai.seg.column
#' @param sample.name
#' @param chr.name
#' @param CN.column
#' @param CN.baseline CN value with amplifications plotted above and deletions plotted below
#' @param CN.track.title
#' @param close.gaps
#' @param close.NA.segs
#' @param CRISPR.column
#' @param RNAi.column
#' @param CRISPR.track.title
#' @param RNAi.track.title
#' @param highlight.gene.col
#' @param highlight.gene.cex
#' @param CRISPR.ylim
#' @param RNAi.ylim
#' @param CN.ylim
#' @param ...
#'
#' @return list of trackList and other Gviz options passed to \code{do.call(plotTracks, .))}
#' @importFrom Gviz DataTrack OverlayTrack IdeogramTrack GenomeAxisTrack
#' @export
#'
cn_dep_gviz <- function(plot.region,
                        crispr, rnai=NULL, cn,
                        highlight.gene=NA,
                        crispr.seg.column=NULL, rnai.seg.column=NULL,
                        sample.name = "",
                        chr.name=NULL,
                        CN.column="CN",
                        CN.baseline=2,
                        CN.track.title="CN",
                        close.gaps=F,
                        close.NA.segs=F,
                        CRISPR.column="Dep",
                        RNAi.column="Dep",
                        CRISPR.track.title="Viability",
                        RNAi.track.title="Viability",
                        CRISPR.seg.col="#984ea3",
                        RNAi.seg.col="#984ea3",
                        CRISPR.seg.lwd=2.5,
                        RNAi.seg.lwd=2.5,
                        highlight.gene.col="#ff7f00",
                        highlight.gene.cex=1.5,
                        CRISPR.ylim=NULL,
                        RNAi.ylim=NULL,
                        CN.ylim=NULL,
                        add_gridlines=T,
                        ...) {

    library(Gviz)

    chr <- as.character(seqnames(plot.region))
    start_pos <- max(1, start(plot.region) - width(plot.region)*.02)
    end_pos <- min(seqlengths(plot.region)[chr], end(plot.region) + width(plot.region)*.02)

    cn$CN <- mcols(cn)[,CN.column]
    if (!is.null(crispr.seg.column)) {
        cn$Seg.Dep <-  mcols(cn)[,crispr.seg.column]
    } else {
        cn$Seg.Dep <- NA
    }

    if (!is.null(rnai.seg.column)) {
        cn$Seg.Dep.RNAi <-  mcols(cn)[,rnai.seg.column]
    } else {
        cn$Seg.Dep.RNAi <- NA
    }

    cn <- cn[seqnames(cn) == chr]


    if (start(cn)[1] != 1) {
        cn <- c(cn[1], cn)
        ranges(cn[1]) <- IRanges(1, start(cn[2])-1)
        cn$CN[1] <- NA
        cn$Seg.Dep[1] <- NA
    }
    if (end(cn)[length(cn)] != seqlengths(cn)[chr]) {
        cn <- c(cn, cn[length(cn)])
        ranges(cn[length(cn)]) <- IRanges(end(cn[length(cn)-1]+1), seqlengths(cn)[chr])
        cn$CN[length(cn)] <- NA
        cn$Seg.Dep[length(cn)] <- NA
    }

    if (close.gaps) {
        gaps <- start(cn)[-1] - end(cn)[-length(cn)] - 1
        end(cn)[-length(cn)] <- end(cn)[-length(cn)] + ceiling(gaps/2)
        start(cn)[-1] <- start(cn)[-1] - floor(gaps/2)
    }

    cn <- cn[cn %over% plot.region]
    start(cn[1]) <- start(plot.region)
    end(cn[length(cn)]) <- end(plot.region)

    cn.starts <- cn

    end(cn.starts) <- start(cn.starts)

    cn.ends <- cn
    start(cn.ends) <- end(cn.ends)

    cn.plot <- c(cn.starts, cn.ends) %>% sort()

    crispr <- crispr[crispr %over% plot.region]

    # crispr
    crispr$CRISPR <- mcols(crispr)[[CRISPR.column]]
    names(crispr) <- crispr$Gene
    crispr$RNAi <- NA_real_

    if (!is.null(rnai)) {
        rnai <- rnai[rnai %over% plot.region]
        rnai$RNAi <- mcols(rnai)[[RNAi.column]]

        names(rnai) <- rnai$Gene
        rnai <- rnai[names(rnai) %in% names(crispr)]
        crispr[rnai$Gene]$RNAi <- mcols(rnai)[[RNAi.column]]
    }

    if (add_gridlines) {
        cn_type <- c("polygon", "g")
        dep_type <- c("p", "g")
    } else {
        cn_type <- "polygon"
        dep_type <- "p"
    }

    if (is.null(CN.ylim)) CN.ylim <- c(min(0, cn.plot$CN, na.rm=T), max(cn.plot$CN, na.rm=T))
    cntrack <- DataTrack(cn.plot, data=cn.plot$CN,
                         legend=F,
                         type=cn_type,
                         ylim=CN.ylim,
                         name=CN.track.title,
                         fill.mountain=c("#377eb8", "#e41a1c"),
                         lty=0,
                         #lwd.mountain=0,
                         baseline=CN.baseline,
                         col.baseline="black",
                         lwd.baseline=2,
                         lty.baseline=1,
                         col.grid="lightgrey", v=-20,
                         col.title="black", col.axis="black")


    if (is.null(CRISPR.ylim)) CRISPR.ylim <- c(min(crispr$CRISPR, na.rm=T), max(0, crispr$CRISPR, na.rm=T))

    crisprtrack <- DataTrack(crispr, data=crispr$CRISPR,
                             name=CRISPR.track.title,
                             type = dep_type,
                             ylim=CRISPR.ylim,
                             cex=0.7, alpha=1, alpha.title=1,
                             legend = TRUE, cex.legend=1,
                             col="black", lwd=2,
                             col.grid="lightgrey", v=-20,
                             col.title="black", col.axis="black")
    overlaytracklist <- list(crisprtrack)

    if (!is.na(highlight.gene)) {
        gene.highlight.crispr <- crispr[crispr$Gene==highlight.gene]
        genehighlighttrack <- DataTrack(gene.highlight.crispr,
                                        data=gene.highlight.crispr$CRISPR,
                                        type = c("p"),
                                        ylim=CRISPR.ylim,
                                        cex=highlight.gene.cex, alpha=1,
                                        alpha.title=1,
                                        legend = TRUE, cex.legend=1,
                                        col=highlight.gene.col, lwd=2,
                                        col.grid="lightgrey", v=-20,
                                        cex.title=1,
                                        col.title="black", col.axis="black")

        overlaytracklist <- c(overlaytracklist, genehighlighttrack)
    }


    if (!is.null(crispr.seg.column)) {

        if (close.NA.segs) {
            cn.plot$Seg.Dep <- fillNAsWithAdjacent(cn.plot$Seg.Dep)
        }

        segtrack <- DataTrack(cn.plot, data=cn.plot$Seg.Dep,
                              type="l", legend=F, lwd=CRISPR.seg.lwd,
                              col=CRISPR.seg.col, ylim=CRISPR.ylim)

        overlaytracklist <- c(overlaytracklist, segtrack)

    }

    crisprtrack <- OverlayTrack(trackList=overlaytracklist,
                                name=CRISPR.track.title)



    if (!is.null(rnai)) {

        if (is.null(RNAi.ylim)) RNAi.ylim <- range(rnai$RNAi, na.rm=T)

        rnaitrack <- DataTrack(rnai, data=rnai$RNAi,
                               name=RNAi.track.title,
                               type = dep_type,
                               ylim=RNAi.ylim,
                               cex=0.7, alpha=1, alpha.title=1,
                               legend = TRUE, cex.legend=1,
                               col="black", lwd=2,
                               col.grid="lightgrey", v=-20,
                               col.title="black", col.axis="black")

        overlaytracklist <- list(rnaitrack)

        if (!is.na(highlight.gene)) {
            # crispr$HighlightGene <- ifelse(crispr$Gene==highlight.gene, crispr$CRISPR, 100)
            gene.highlight.rnai <- rnai[rnai$Gene==highlight.gene]
            genehighlighttrack <- DataTrack(gene.highlight.rnai,
                                            data=gene.highlight.rnai$RNAi,
                                            type = c("p"),
                                            ylim=RNAi.ylim,
                                            cex=highlight.gene.cex, alpha=1,
                                            alpha.title=1,
                                            legend = TRUE, cex.legend=1,
                                            col=highlight.gene.col, lwd=2,
                                            col.grid="lightgrey", v=-20,
                                            cex.title=1,
                                            col.title="black", col.axis="black")

            overlaytracklist <- c(overlaytracklist, genehighlighttrack)
        }


        if (!is.null(rnai.seg.column)) {

            if (close.NA.segs) {
                cn.plot$Seg.Dep.RNAi <- fillNAsWithAdjacent(cn.plot$Seg.Dep.RNAi)
            }

            segtrack <- DataTrack(cn.plot, data=cn.plot$Seg.Dep.RNAi,
                                  type="l", legend=F, lwd=RNAi.seg.lwd,
                                  col=RNAi.seg.col, ylim=RNAi.ylim)

            overlaytracklist <- c(overlaytracklist, segtrack)

        }

        rnaitrack <- OverlayTrack(trackList=overlaytracklist, name=RNAi.track.title)
    }

    ideotrack <- IdeogramTrack(chromosome=chr,
                               genome="hg19",
                               name=chr.name,
                               showId=TRUE, cex=1,
                               fontsize=12,
                               fontcolor="black",
                               fontfamily="sans",
                               showBandId = FALSE)
    blanktrack <- GenomeAxisTrack(col="transparent", fontcolor="transparent")
    axistrack <- GenomeAxisTrack(labelPos="above",
                                 col="black", fontcolor="black")



    if (is.null(rnai)) {
        return(list(trackList=list(ideotrack, blanktrack, axistrack, cntrack, blanktrack, crisprtrack),
                    sizes=c(1, 0.5, 1.5, 4, 0.5, 6),
                    from=start_pos, to=end_pos,
                    main = sample.name,
                    background.title="transparent",
                    col.border.title="transparent", ...))
    } else {
        return(list(trackList=list(ideotrack, blanktrack, axistrack, cntrack, blanktrack,
                                   crisprtrack, blanktrack, rnaitrack),
                    sizes=c(1, .5, 1.5, 4, 0.5, 6, 0.5, 6),
                    from=start_pos, to=end_pos,
                    main = sample.name,
                    background.title="transparent",
                    col.border.title="transparent", ...))
    }

}

fillGaps <- function(gr) {
    gr.list <- split(gr, seqnames(gr))
    llply(gr.list, function(g) {
        gaps <- start(g)[-1] - end(g)[-length(g)] - 1
        end(g)[-length(g)] <- end(g)[-length(g)] + ceiling(gaps/2)
        start(g)[-1] <- start(g)[-1] - floor(gaps/2)
        return(g)
    }) %>% GRangesList() %>% unlist() %>% set_names(names(gr))
}

fillNAsWithAdjacent <- function(x) {

    if (length(x) < 3) return(x)
    if (!any(is.na(x))) return(x)

    x.na <- is.na(x)
    x.na.rle <- rle(x.na)
    x.na.rle$values[1] <- NA
    x.na.rle$values[length(x.na.rle$values)] <- NA

    replacement.rle <- rep(NA, times=length(x.na.rle$values))

    left.ind <- cumsum(x.na.rle$lengths) - x.na.rle$lengths
    right.ind <- cumsum(x.na.rle$lengths) + 1

    left.ind <- left.ind[2:(length(left.ind)-1)]
    right.ind <- right.ind[2:(length(right.ind)-1)]


    replacement.rle[2:(length(replacement.rle)-1)] <- ifelse(abs(x[left.ind]) <
                                                                 abs(x[right.ind]),
                                                             x[left.ind], x[right.ind])
    x.na.rle$values <- replacement.rle

    replacement.vec <- inverse.rle(x.na.rle)

    x[!is.na(replacement.vec)] <- replacement.vec[!is.na(replacement.vec)]

    return(x)

}
