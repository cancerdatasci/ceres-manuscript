
#' Get a PAM sequence
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome getSeq
#' @export
#'
fetchPAM <- function(chr, pos, strand) {
    PAM_start <- ifelse(strand=="+", pos+20, pos-3)
    PAM_end <- ifelse(strand=="+", pos+22, pos-1)

    PAM <- getSeq(BSgenome.Hsapiens.UCSC.hg19,
                  names=chr,
                  start=PAM_start,
                  end=PAM_end,
                  strand=strand,
                  as.character=TRUE)
    return(PAM)
}




#' Read BAM file report, filter with valid PAM sequence, report # of alignments
#' @param bam.file BAM file path
#' @param max.alns only consider guides with fewer than this many alignments
#' @param pam PAM sequence as regular expression
#' @param include.no.align include guides that do not map to genome
#' @param as.df return data.frame if \code{TRUE}, GRanges object if \code{FALSE}
#' @return data.frame or GRanges object of results
#' @import dplyr
#' @importFrom Rsamtools scanBam
#' @importFrom GenomeInfoDb Seqinfo
#' @importFrom BSgenome.Hsapiens.UCSC.hg19 BSgenome.Hsapiens.UCSC.hg19
#' @importFrom BSgenome getSeq
#' @export
#'
guideAlignments <- function(bam.file, max.alns=100, pam="[ACGTN]GG",
                            include.no.align=F, as.df=T) {
    chromosomes <- c(as.character(1:22), "X", "Y") %>% str_c("chr", .)

    bam <- Rsamtools::scanBam(bam.file)

    ### Bam file to data frame, count number of alignments
    bamdf <- as_data_frame(bam[[1]][1:6]) %>%
        mutate(Gene = str_extract(qname, "_.+$") %>% str_sub(start=2),
               Guide = str_extract(qname, "^[A-Z]+"),
               Aligned = !is.na(rname)) %>%
        filter(Aligned & rname %in% chromosomes) %>%
        group_by(qname) %>%
        mutate(N.alns = sum(Aligned)) %>%
        ungroup()

    bamdf.noAlign <- as_data_frame(bam[[1]][1:6]) %>%
        mutate(Gene = str_extract(qname, "_.+$") %>% str_sub(start=2),
               Guide = str_extract(qname, "^[A-Z]+"),
               Aligned = !is.na(rname)) %>%
        filter(!Aligned) %>%
        mutate(N.alns = 0)


    ### Limit guides with too many alignemnts, annotate position from biomart data
    bamdf.filt <- bamdf %>% filter(N.alns < max.alns)

    bamdf.pam <- mutate(bamdf.filt,
                        Cut.Pos = ifelse(strand=="+", pos+16, pos+2),
                        PAM.start = ifelse(strand=="+", pos+20, pos-3),
                        PAM.end = ifelse(strand=="+", pos+22, pos-1))

    bamdf.pam$PAM <-
        BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
                         names=bamdf.pam$rname,
                         start=bamdf.pam$PAM.start,
                         end=bamdf.pam$PAM.end,
                         strand=bamdf.pam$strand,
                         as.character=TRUE)

    bamdf.pam <- bamdf.pam %>%
        group_by(qname) %>%
        mutate(N.alns = sum(str_detect(PAM, pam))) %>%
        ungroup()


    if (any(bamdf.pam$N.alns == 0)) {
        bamdf.noAlign <- bind_rows(bamdf.noAlign,
                                   bamdf.pam %>%
                                       filter(N.alns == 0) %>%
                                       distinct(qname, Gene, Guide, .keep_all=T) %>%
                                       select(qname, Gene, Guide, N.alns) %>%
                                       mutate(Aligned = FALSE))
        bamdf.pam <- bamdf.pam %>% filter(str_detect(PAM, pam))
    }


    if (include.no.align) {
        bamdf.final <- bind_rows(bamdf.pam, bamdf.noAlign)
    } else {
        bamdf.final <- bamdf.pam
    }

    if (as.df) {
        return(bamdf.final)
    } else {
        hg19info <- Seqinfo(genome="hg19")[chromosomes]

        return(makeGRangesFromDataFrame(bamdf.final,
                                        keep.extra.columns=T,
                                        seqinfo=hg19info))
    }
}
