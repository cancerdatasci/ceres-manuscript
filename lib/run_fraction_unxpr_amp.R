

fraction_amps_unxpr <- function(gene_level_df, amp_cutoff=4, rpkm_cutoff=-1,
                                xlims=c(0, -8)) {

    fp_df <- gene_level_df %>%
        dplyr::rename(CERES = Gene_Effect,
                      `Uncorrected Data` = AvgGuide) %>%
        mutate(Amp = CN > amp_cutoff,
               XPR = RPKM > rpkm_cutoff) %>%
        filter(Gene !=  "", !is.na(Amp), !is.na(XPR), !is.na(`Uncorrected Data`), !is.na(CERES)) %>%
        gather(Method, Dep, `Uncorrected Data`, CERES) %>%
        mutate(Method = factor(Method, levels=c("Uncorrected Data", "CERES"), ordered=T)) %>%
        group_by(Method, Gene) %>%
        mutate(RelDep = Dep - mean(Dep, na.rm=T)) %>%
        ungroup() %>%
        group_by(Method, CellLine) %>%
        mutate(RelDep = (RelDep - mean(RelDep)) / sd(RelDep)) %>%
        ungroup %>%
        group_by(Method, Gene) %>%
        filter(any(Dep < 0)) %>%
        ungroup



    fp_amps_rel <- fp_df %>%
        arrange(Method, RelDep) %>%
        group_by(Method) %>%
        mutate(FracDepAmp = cumsum(Amp)/(1:length(Amp)),
               ScaledRank = percent_rank(RelDep)) %>%
        ungroup

    fp_xpr_rel <- fp_df %>%
        arrange(Method, RelDep) %>%
        group_by(Method) %>%
        mutate(FracDepUnXPR = cumsum(!XPR)/(1:length(XPR)),
               ScaledRank = percent_rank(RelDep)) %>%
        ungroup




    fp_amps_rel_sample <- sample_fp_df(fp_amps_rel, data_column=RelDep)

    fp_xpr_rel_sample <- sample_fp_df(fp_xpr_rel, data_column=RelDep)



    frac_amps_gg <- ggplot(fp_amps_rel_sample, aes(x=RelDep, y=FracDepAmp, color=Method)) +
        geom_line(size=1) +
        scale_x_reverse(limits=xlims, breaks=seq(0,-10, by=-2), expand=c(0,0)) +
        scale_y_continuous(limits=c(0, 0.4), breaks=seq(0,0.5,by=0.1),
                           expand=c(0,0), labels=scales::percent) +
        scale_color_manual(values=c("Uncorrected Data" = "#e41a1c","CERES" = "#377eb8")) +
        labs(y="Amplified Genes (%)",
             x="Differential Dependency (Z-Score)") +
        theme(legend.position=c(0.05, 0.95),
              legend.justification=c(0,1),
              legend.title=element_blank())

    frac_unxpr_gg <- ggplot(fp_xpr_rel_sample, aes(x=RelDep, y=FracDepUnXPR, color=Method)) +
        geom_line(size=1) +
        scale_x_reverse(limits=xlims, breaks=seq(0, -10, by=-2), expand=c(0,0)) +
        scale_y_continuous(limits=c(0, .4), breaks=seq(0,0.5,by=0.1),
                           expand=c(0,0), labels=scales::percent) +
        scale_color_manual(values=c("Uncorrected Data" = "#e41a1c","CERES" = "#377eb8")) +
        labs(y="Unexpressed Genes (%)",
             x="Differential Dependency (Z-Score)") +
        theme(legend.position=c(0.95, 0.95),
              legend.justification=c(1,1),
              legend.title=element_blank())

    return(list(frac_amps_df=fp_amps_rel, frac_xpr_df=fp_xpr_rel,
                frac_amps_gg=frac_amps_gg, frac_unxpr_gg=frac_unxpr_gg))
}



sample_fp_df <- function(df, data_column=Dep) {

    data_column <- enquo(data_column)

    sample_df <- df %>%
        mutate(RANK = min_rank(!!data_column)) %>%
        {bind_rows(filter(., RANK < 1e3),
                   filter(., RANK >= 1e3, RANK < 1e4) %>% group_by(Method) %>%
                       sample_n(1000, replace=T) %>% ungroup(),
                   filter(., RANK >= 1e4, RANK < 1e5) %>% group_by(Method) %>%
                       sample_n(500, replace=T) %>% ungroup(),
                   filter(., RANK >= 1e5, RANK < 1e6) %>% group_by(Method) %>%
                       sample_n(200, replace=T) %>% ungroup(),
                   filter(., RANK >= 1e6) %>% group_by(Method) %>%
                       sample_n(100, replace=T) %>% ungroup())} %>%
        select(-RANK)

    return(sample_df)
}
