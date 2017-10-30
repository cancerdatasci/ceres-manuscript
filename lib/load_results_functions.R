
load_ceres_outputs <- function(output_dir) {

    if (!dir.exists(output_dir)) {
        warning("Warning: output directory does not exist")
        return(NULL)
    }

    gene_effect_fit <-
        readRDS(file.path(output_dir, "ge_fit.rds"))

    cutting_effect_fit <-
        readRDS(file.path(output_dir, "ce_fit.rds"))

    guide_activity_fit <-
        readRDS(file.path(output_dir, "sgrna_fit.rds"))

    fit_summary <-
        read.table(file.path(output_dir, "summary.txt"), header=T)

    ceres_output <- list(gene_effect = gene_effect_fit %>%
                             mat.to.df("Gene", "CellLine", "Gene_Effect"),
                         guide_activity = guide_activity_fit %>%
                             select(Guide = sgRNA,
                                    Activity=efficacy,
                                    Offset = offset),
                         cutting_effect = cutting_effect_fit %>%
                             filter(position == "end") %>%
                             select(CellLine = cell_line,
                                    segment, slope, intercept),
                         summary = fit_summary)
}

load_ceres_results <- function(ceres_run_date, sub_dir, input_data,
                               base_dir="~/Google Drive/CERES/ModelOutputs",
                               new_ceres=T) {

    fit_dir <- file.path(base_dir, ceres_run_date, sub_dir)
    if (!dir.exists(fit_dir)) return(NULL)


    if (!file.exists(file.path(fit_dir, "sgRNA.effect.fitted.rds"))) {
        sgRNA_effect_fitted <- NULL
    } else {
        sgRNA_effect_fitted <-
            readRDS(file.path(fit_dir, "sgRNA.effect.fitted.rds"))
    }
    cutting_intercept_fitted <-
        readRDS(file.path(fit_dir, "intercept.fitted.rds"))
    cutting_toxicity_fitted <-
        readRDS(file.path(fit_dir, "cutting.toxicity.fitted.rds"))
    gene_effect_fitted <-
        readRDS(file.path(fit_dir, "gene.effect.fitted.rds"))
    hinge_point_fitted <-
        readRDS(file.path(fit_dir, "hinge.point.fitted.rds"))
    guide_efficacy_fitted <-
        readRDS(file.path(fit_dir, "guide.prob.ko.locus.fitted.rds")) %>%
        as_data_frame()

    if (file.exists(file.path(fit_dir, "quantiles.rds"))) {
        quantiles <- readRDS(file.path(fit_dir, "quantiles.rds"))
    } else {
        quantiles <- NULL
    }

    if (file.exists(file.path(fit_dir, "guide.offset.fitted.rds"))) {
        guide_offset_fitted <-
            readRDS(file.path(fit_dir, "guide.offset.fitted.rds"))
    } else {
        guide_offset_fitted <- NULL
    }

    # validation_set <- readRDS(file.path(fit_dir,
    #                                     str_c("validationSet.rds")))


    guide_efficacy <- guide_efficacy_fitted %>%
        select(Guide, Locus, Efficacy=Prob)

    # if (is.null(guide_offset_fitted)) guide_offset_fitted <- 0

    # if (is.matrix(guide_offset_fitted) && ncol(guide_offset_fitted) == 3) {
    #     guide_offset <- data_frame(Guide=guide_efficacy$Guide,
    #                                Offset=guide_offset_fitted[,3]) %>%
    #         distinct()
    # } else if (is.vector(guide_offset_fitted)) {
    #     guide_offset <- data_frame(Guide=unique(c(input_data$sgRNA_effect_observed$Guide,
    #                                               input_data$negctrl_sgRNA_observed$Guide)),
    #                                Offset=guide_offset_fitted)
    # } else {
    #     stop("don't know how to handle the guide offsets\n")
    # }

    if (new_ceres) {
        sgRNA_fitted_df <- NULL
    } else {
        sgRNA_fitted_df <- sgRNA_effect_fitted %>%
            as.data.frame() %>%
            mutate(Guide = rownames(.)) %>%
            gather(CellLine, sgRNA_Fitted, -Guide)
    }

    if (is.null(cutting_intercept_fitted)) cutting_intercept_fitted <- 0

    if (new_ceres) {

        cutting_intercept_df <- data_frame(CellLine = colnames(gene_effect_fitted),
                                           Intercept = cutting_intercept_fitted)
        cutting_slopes_df <- cutting_toxicity_fitted %>%
            set_rownames(colnames(gene_effect_fitted)) %>%
            set_colnames(1:ncol(.)) %>%
            as.data.frame() %>%
            mutate(CellLine = rownames(.)) %>%
            gather(Segment, Slope, -CellLine)

        cutting_starts_df <-quantiles[,seq(1, ncol(quantiles)-1, by=2)] %>%
            set_rownames(colnames(gene_effect_fitted)) %>%
            set_colnames(1:ncol(.)) %>%
            as.data.frame() %>%
            mutate(CellLine = rownames(.)) %>%
            gather(Segment, Start, -CellLine)
        cutting_hinges_df <- quantiles[,seq(2, ncol(quantiles), by=2)] %>%
            set_rownames(colnames(gene_effect_fitted)) %>%
            set_colnames(1:ncol(.)) %>%
            as.data.frame() %>%
            mutate(CellLine = rownames(.)) %>%
            gather(Segment, Hinge, -CellLine)
        segment_intercepts_df <-
            (cutting_toxicity_fitted *
            (quantiles[,seq(2, ncol(quantiles), by=2)] - quantiles[,seq(1, ncol(quantiles)-1, by=2)])) %>%
            aaply(1, function(x) cumsum(c(0,head(x,n=-1)))) %>%
            set_rownames(colnames(gene_effect_fitted)) %>%
            set_colnames(1:ncol(.)) %>%
            as.data.frame() %>%
            mutate(CellLine = rownames(.)) %>%
            gather(Segment, SegIntercept, -CellLine)

        cutting_sensitivity_df <- cutting_intercept_df %>%
            left_join(cutting_starts_df) %>%
            left_join(cutting_hinges_df) %>%
            left_join(segment_intercepts_df) %>%
            left_join(cutting_slopes_df) %>%
            mutate(Segment = as.numeric(Segment))


    } else {
        cutting_sensitivity_df <- data_frame(CellLine = colnames(gene_effect_fitted),
                                             CERES_Intercept = cutting_intercept_fitted[1:nrow(cutting_toxicity_fitted)],
                                             CERES_Slope = cutting_toxicity_fitted[,1],
                                             CERES_Hinge = hinge_point_fitted)# / cutting_toxicity_fitted[,1])
    }


    if (is.null(colnames(gene_effect_fitted))) {
        colnames(gene_effect_fitted) <- colnames(sgRNA_effect_fitted)
    }

    if (is.null(rownames(gene_effect_fitted))) {
        rownames(gene_effect_fitted) <- unique(input_data$locus_gene$Gene)
    }


    gene_effect_df <- gene_effect_fitted %>%
        as.data.frame() %>%
        mutate(Gene = rownames(.)) %>%
        gather(CellLine, Gene_Effect, -Gene)


    ceres_output <- list(gene_effect = gene_effect_df,
                         sgRNA_effect_fitted = sgRNA_fitted_df,
                         guide_efficacy = guide_efficacy,
                         # guide_offset = guide_offset,
                         cutting_sensitivity = cutting_sensitivity_df)
}


load_ceres_log <- function(ceres_run_date, fit_suffix,
                           base_dir="~/Google\ Drive/CERES/ModelOutputs"){

    fit_dir <- file.path(base_dir, ceres_run_date, fit_suffix, "log")
    if (!dir.exists(fit_dir)) return(NULL)

    log_tab <- readr::read_tsv(file.path(fit_dir, str_c("logfile_", fit_suffix, ".txt")))

}

load_ceres_optimization <- function(ceres_run_date, fit_suffix,
                                    base_dir="~/Google\ Drive/CERES/ModelOutputs"){

    fit_dir <- file.path(base_dir, ceres_run_date, fit_suffix, "optimization")
    if (!dir.exists(fit_dir)) return(NULL)

    errors <- c()
    lg <- c()
    lq <- c()
    train_errors <- c()

    for(file in dir(file.path(fit_dir))[grepl("job_", dir(file.path(fit_dir)))]){

        filepath <- file.path(fit_dir, file)

        errors <- c(errors, read_tsv(filepath) %>% magrittr::extract2(., "test_error"))
        lg <- c(lg, read_tsv(filepath) %>% magrittr::extract2(., "lg"))
        lq <- c(lq, read_tsv(filepath) %>% magrittr::extract2(., "lq"))

    }

    for(file in dir(file.path(fit_dir, "log"))[grepl(".txt", dir(file.path(fit_dir, "log")))]){
        filepath <- file.path(fit_dir, "log", file)
        train_errors <- c(train_errors, read_tsv(filepath) %>% magrittr::extract2(., "train_err") %>% min())
    }

    opt_tab <- data.frame(lg=lg, lq=lq, test_error=errors, train_error=train_errors)

}


load_ceres_inputs <- function(inputs_dir,
                              replicates=F, do_parallel=F, collapse_replicates=T) {

    # dir.create(inputs_dir, recursive=T, showWarnings=F)


    # if (!is.null(prepare_inputs_fun)) prepare_inputs_fun(inputs_dir, ...)

    guide_dep <- readRDS(file.path(inputs_dir, "guide_sample_dep.Rds"))
    if (file.exists(file.path(inputs_dir, "negctrl_guide_sample_dep.Rds"))) {
        negctrl_dep <- readRDS(file.path(inputs_dir, "negctrl_guide_sample_dep.Rds"))
    } else {
        negctrl_dep <- NULL
    }
    guide_locus_df <- readRDS(file.path(inputs_dir, "guide_locus_map3.Rds"))
    locus_gene_df <- readRDS(file.path(inputs_dir, "locus_gene_map3.Rds"))
    locus_cn <- readRDS(file.path(inputs_dir, "locus_sample_cn.Rds"))
    gene_ids <- readRDS(file.path(inputs_dir, "gene_id_table.Rds"))

    if (replicates) {
        replicate_map <- readRDS(file.path(inputs_dir, "replicate_map.Rds")) %>%
            select(Replicate, CellLine)
    } else {
        replicate_map <- NULL
    }

    locus_cn_df <- locus_cn %>%
        as.data.frame() %>%
        mutate(Locus = rownames(locus_cn)) %>%
        gather(CellLine, CN, -Locus)

    if (replicates && collapse_replicates) {
        sgRNA_observed_df <-
            collapse_reps(guide_dep, replicate_map, do_parallel=do_parallel) %>%
            as.data.frame() %>%
            mutate(Guide = rownames(.)) %>%
            gather(CellLine, sgRNA_Observed, -Guide)

        if (!is.null(negctrl_dep)) {
            negctrl_sgRNA_observed_df <-
                collapse_reps(negctrl_dep, replicate_map, do_parallel=do_parallel) %>%
                as.data.frame() %>%
                mutate(Guide = rownames(.)) %>%
                gather(CellLine, sgRNA_Observed, -Guide)
        } else {
            negctrl_sgRNA_observed_df <- NULL
        }

    } else {
        sgRNA_observed_df <- guide_dep %>%
            as.data.frame() %>%
            mutate(Guide = rownames(.)) %>%
            gather(CellLine, sgRNA_Observed, -Guide)

        if (!is.null(negctrl_dep)) {
            negctrl_sgRNA_observed_df <- negctrl_dep %>%
                as.data.frame() %>%
                mutate(Guide = rownames(.)) %>%
                gather(CellLine, sgRNA_Observed, -Guide)
        } else {
            negctrl_sgRNA_observed_df <- NULL
        }
    }


    ceres_input <- list(
        sgRNA_effect_observed = sgRNA_observed_df,
        negctrl_sgRNA_observed = negctrl_sgRNA_observed_df,
        guide_locus = guide_locus_df %>% select(Guide, Locus) %>% distinct(),
        locus_gene = locus_gene_df %>% select(Locus, Gene) %>% distinct(),
        locus_cn = locus_cn_df,
        gene_ids = gene_ids %>% distinct(GeneID, .keep_all=T),
        replicate_map = replicate_map)

}

