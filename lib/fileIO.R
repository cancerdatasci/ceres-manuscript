#' read a *.gct file into your session
#' @param gct.file path to gct formatted file
#' @return matrix with a "description" attribute
#' @importFrom magrittr set_rownames
#' @export
read.gct <- function(gct.file) {
    dat <- read.delim(gct.file, skip=2, header=TRUE, sep="\t",
                      check.names=F, as.is=T)
    if (colnames(dat)[1] == "NAME") {
        dat %<>% rename(Name = NAME)
    }
    if (colnames(dat)[2] == "DESCRIPTION") {
        dat %<>% rename(Description = DESCRIPTION)
    }
    mat <- as.matrix(select(dat, -Name, -Description)) %>%
        set_rownames(dat$Name)
    attr(mat,"Description") <- dat$Description %>% set_names(dat$Name)
    return(mat)
}


#' write a *.gct file to disk
#' @param mat data matrix, optionally with attribute "description"
#' @param gct.file path to file to write
#' @param description character vector; length must equal number of rows in mat
#' @export
write.gct <- function(mat, gct.file, description = attr(mat, "Description")) {

    stopifnot(!is.null(rownames(mat)))
    stopifnot(!is.null(colnames(mat)))

    if(is.null(description)) {
        description <- rownames(mat)
    }

    stopifnot(length(description) == nrow(mat))


    dat <- as.data.frame(mat) %>%
        cbind(Name=rownames(mat), Description=description, .)

    cat("#1.2\n", file=gct.file)
    cat(nrow(mat), ncol(mat), sep="\t", file=gct.file, append = TRUE)
    cat("\n", file=gct.file, append =  TRUE)
    cat(colnames(dat), sep="\t", file=gct.file, append = TRUE)
    cat("\n", file=gct.file, append =  TRUE)
    suppressWarnings(
        write.table(dat, file=gct.file, append = TRUE,
                    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE )
    )

}

#' read a *.gmt file into your session
#' @param gmt.file path to gmt formatted file
#' @param as.df whether to return as data frame or list
#' @return either a list of character vectors or a two column data.frame
#' @export
read.gmt <- function(gmt, as.df=T) {
    genesets <- readLines(gmt, warn=F) %>%
        str_split("\t") %>%
        set_names(., laply(., function(l) l[1])) %>%
        llply(function(l) l[3:length(l)])

    if (as.df) {
        genesets.df <- ldply(genesets,
                             function (gs) data_frame(Gene = gs),
                             .id="GeneSet")
        return(genesets.df)
    } else {
        return(genesets)
    }
}
