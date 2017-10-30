

run_uncorrected_slopes <- function(ceres_input,
                                   neg_ctrl_guides, pos_ctrl_guides,
                                   intercept=T,
                                   do_parallel,
                                   draw_plots=T, plot_dir=".",
                                   test_breakpoints=c(seq(0, 10, by=2), seq(20, 50, by=10)),
                                   cut_bins=seq(-1, 51, by=2),
                                   dens_binwidth_x=0.25, dens_binwidth_y=0.05,
                                   dep_ylims=c(-6,2)) {


    cat("Tidying data...\n")

    guide_cuts <-
        collapse_rows_of_matrix(mat=ceres_input$locus_cn %>%
                                    select(Locus, CellLine, CN) %>%
                                    df.to.mat(),
                                group_df=ceres_input$guide_locus,
                                collapse_fun=Matrix::colSums, na.rm=T,
                                group_var=Guide, sample_var=Locus, do_parallel=do_parallel) %>%
        rbind(matrix(0, nrow=length(neg_ctrl_guides$Guide),  ncol=ncol(.),
                     dimnames=list(neg_ctrl_guides$Guide, colnames(.))))

    uncorrected_dep_cuts <-
        make.a.tidy.dataset(list(Dep=ceres_input$sgRNA_effect_observed %>%
                                     select(Guide, CellLine, sgRNA_Observed) %>%
                                     df.to.mat(),
                                 Cuts=guide_cuts),
                            dim.names=c("Guide", "CellLine")) %>%
        filter(!is.na(Dep), !is.na(Cuts)) %>% distinct


    cat("Running test breakpoints...\n")

    cuts_pwlm_rsquareds <- uncorrected_dep_cuts %>%
        ddply(.(CellLine), function(d) {
            lms <- test_saturating_lm(d$Cuts, d$Dep, intercept=intercept,
                                      testpoints=test_breakpoints, do.parallel=F)
            df <- data_frame(breakpoint=as.numeric(names(lms)),
                             LM = lms)
            df0 <- df %>% filter(breakpoint == 0) %>%
                mutate(r.squared = laply(LM, function(x) summary(x)$r.squared),
                       adj.r.squared = r.squared) %>%
                select(-LM)
            df1 <- df %>% filter(breakpoint > 0) %>% rowwise() %>% glance(LM)
            bind_rows(df0, df1)
        }, .parallel=do_parallel) %>% as_data_frame()

    best_pwlm <- cuts_pwlm_rsquareds %>% group_by(CellLine) %>%
        summarise(breakpoint = breakpoint[r.squared == max(r.squared)][1],
                  r.squared = max(r.squared)) %>%
        ungroup

    cat("Optimizing breakpoint models...\n")
    cuts_pwlm <- uncorrected_dep_cuts %>%
        ddply(.(CellLine), function(d) {
            initial_breakpoint <- best_pwlm %>%
                filter(CellLine==d$CellLine[1]) %$%
                breakpoint
            fit_saturating_lm(d$Cuts, d$Dep, intercept=intercept, init.par=initial_breakpoint)
        }, .parallel=do_parallel) %>% as_data_frame()


    fit_pwlm_df <- cuts_pwlm %>% rowwise() %>%
        augment(model, newdata=data_frame(x=0:50)) %>%
        ungroup %>%
        select(CellLine, x, .fitted, .se.fit)

    cuts_pwlm <- cuts_pwlm %>% select(-model)


    cut_bins_means_df <- uncorrected_dep_cuts %>%
        filter(Cuts > 0) %>%
        mutate(Cut_Bin = cut(Cuts, cut_bins, labels=seq(0, 50, by=2)) %>%
                   as.character %>% as.numeric()) %>%
        group_by(CellLine, Cut_Bin) %>%
        summarise(Dep = mean(Dep, na.rm=T),
                  N = n()) %>%
        filter(N >= 25, !is.na(Cut_Bin)) %>%
        ungroup



    if (draw_plots) {
        dir.create(plot_dir, showWarnings=F, recursive=T)
        cat("Calculating density...\n")
        uncorrected_dep_cuts <- uncorrected_dep_cuts %>%
            ddply(.(CellLine), function(d) {
                xbreaks <- seq(0, max(d$Cuts)+dens_binwidth_x, by=dens_binwidth_x)
                ybreaks <- seq(min(d$Dep), max(d$Dep)+dens_binwidth_y, by=dens_binwidth_y)
                h2 <- hist2(d$Cuts, d$Dep, xbreaks=xbreaks, ybreaks=ybreaks, plot=FALSE)
                xind <- cut(d$Cuts, breaks=xbreaks, labels=seq(dim(h2$z)[1]), include.lowest=T, right=F)
                yind <- cut(d$Dep, breaks=ybreaks, labels=seq(dim(h2$z)[2]), include.lowest=T, right=F)
                d$Dens <- h2$z[cbind(xind, yind)]
                return(d)
            }, .parallel=do_parallel)
        cat("Drawing figures...\n")
        gg_figs <- dlply(uncorrected_dep_cuts, .(CellLine), function(d) {
            saturating_plot(d, cut_bins_df=cut_bins_means_df,
                            pwlm_df=cuts_pwlm, pwlm_fit=fit_pwlm_df,
                            rsquareds_df=cuts_pwlm_rsquareds,
                            neg_ctrls=neg_ctrl_guides, pos_ctrls=pos_ctrl_guides,
                            dep_ylims=dep_ylims)
        }, .parallel=F)
        l_ply(names(gg_figs), function(cell_line) {
            ggsave(file.path(plot_dir, str_c(cell_line, ".png")),
                   plot=gg_figs[[cell_line]], width=6, height=4)
        }, .parallel=F)
    } else {
        gg_figs=NULL
    }

    return(list(tdat=uncorrected_dep_cuts, pwlm=cuts_pwlm, pwlm_fit=fit_pwlm_df, test_pwlms=cuts_pwlm_rsquareds,
                bin_means=cut_bins_means_df, figs=gg_figs))

}



fit_saturating_lm <- function (x, y, init.par=6, intercept=F) {

    f <- function (breakpoint) {
        hinge <- function(x) ifelse(x < breakpoint, x, breakpoint)
        if (intercept) {
            lm( y ~ hinge(x) )
        } else {
            lm( y ~ hinge(x) + 0 )
        }

    }

    r2 <- function(breakpoint) -1*summary(f(breakpoint))$r.squared

    optim_r2 <- optim(init.par, r2, method="L-BFGS-B", lower=0, upper=max(x))

    best_fit <- f(optim_r2$par)

    return(data_frame(intercept=ifelse(intercept, coef(best_fit)[1], 0),
                      slope=ifelse(intercept, coef(best_fit)[2], coef(best_fit)[1]),
                      breakpoint=optim_r2$par,
                      floor=ifelse(intercept,
                                   coef(best_fit)[1]+coef(best_fit)[2]*optim_r2$par,
                                   coef(best_fit)[1]*optim_r2$par),
                      rsquared=summary(best_fit)$r.squared,
                      convergence=optim_r2$convergence,
                      message=optim_r2$message,
                      model=list(best_fit))
    )

}





test_saturating_lm <- function (x, y, intercept=F, testpoints=NULL,
                                lims=NULL, by=ifelse(is.null(lims), 5, 1),
                                do.parallel=F) {

    f <- function (breakpoint) {
        hinge <- function(x) ifelse(x < breakpoint, x, breakpoint)
        if (intercept) {
            lm( y ~ hinge(x) )
        } else {
            lm( y ~ hinge(x) + 0 )
        }

    }

    r2 <- function(breakpoint) -1*summary(f(breakpoint))$r.squared

    if (is.null(lims)) lims <- c(floor(min(x)), ceiling(max(x)))
    if (is.null(testpoints)) testpoints <- seq(lims[1], lims[2], by=by)

    lms <-  testpoints %>% as.list %>% set_names(.,.) %>%
        llply(f, .parallel=do.parallel)

    return(lms)

}


saturating_plot <- function(d, cut_bins_df, pwlm_df, pwlm_fit, rsquareds_df, neg_ctrls, pos_ctrls,
                            dep_ylims=c(-6,2)) {

    cell_line <- d$CellLine[1]

    cut_bins_df <- cut_bins_df %>% filter(CellLine == cell_line)
    pwlm_df <- pwlm_df %>% filter(CellLine == cell_line)
    pwlm_fit <- pwlm_fit %>% filter(CellLine == cell_line)
    rsquareds_df <- rsquareds_df %>% filter(CellLine == cell_line)

    samp_tdat <- downsample_data(d, "Dens", c(0, 0.1, 0.25, 1), c(1, 0.5, 0.1))

    gg_rsquareds <- rsquareds_df %>% filter(CellLine==cell_line) %>%
        ggplot(aes(x=breakpoint, y=r.squared)) +
        geom_line() +
        geom_vline(aes(xintercept=breakpoint), linetype=2,
                   data=pwlm_df, show.legend=F) +
        theme(axis.title.x = element_text(size=8, margin=margin(1,0,0,0)),
              axis.title.y = element_text(size=8, margin=margin(0,0,2,0)),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_text(size=6, margin=margin(1,0,0,0)),
              plot.background = element_rect(fill=NA, linetype=1, color="black", size=0))


    d_ctrls <- d %>% mutate(Ctrl = case_when(Guide %in% neg_ctrls$Guide ~ "Neg",
                                             Guide %in% pos_ctrls$Guide ~ "Pos",
                                             TRUE ~ "None") %>%
                                factor(levels=c("None", "Neg", "Pos"), ordered=T))
    d_ctrls_means <- d_ctrls %>%
        group_by(Ctrl) %>%
        summarise(Median = median(Dep, na.rm=T),
                  Mean = mean(Dep, na.rm=T))


    gg_scatter <-
        samp_tdat %>%
        arrange(Dens) %>%
        ggplot(aes(x=Cuts, y=Dep, color=Dens)) +
        geom_hline(yintercept=d_ctrls_means$Median, alpha=0.5) +
        geom_point(alpha=0.75, size=1, stroke=0) +
        geom_segment(aes(x=Cut_Bin-0.6, xend=Cut_Bin+0.6, y=Dep, yend=Dep),
                     data=cut_bins_df,
                     color="black", size=1.5) +
        geom_line(aes(x=x, y=.fitted),
                  data=pwlm_fit,
                  color="#e41a1c", alpha=0.75, size=2) +
        scale_color_gradientn("", trans="log10", guide=FALSE,
                              colours = rev(rainbow(10, v=.9, s=0.9, end = 4/6))) +
        coord_cartesian(xlim=c(0, 50), ylim=dep_ylims) +
        labs(x="cuts", y="log2FC",
             title=str_extract(cell_line, "^[^_]+"))

    gg_scatter_gtable <- gg_scatter %>%
        ggplot_build() %>% ggplot_gtable()

    gg_ctrls <- ggplot(d_ctrls %>% filter(!is.na(Dep)),
                       aes(x=Dep, y=..density.., fill=Ctrl)) +
        geom_density(alpha=0.75, show.legend=F) +
        geom_vline(xintercept=d_ctrls_means$Median) +
        scale_x_continuous(limits=dep_ylims) +
        scale_fill_manual(values=c("Pos"="#e41a1c", "Neg"="#377eb8", "None"="grey80")) +
        # labs(x=NULL, y=NULL) +
        coord_flip() +
        theme_minimal() +
        theme(axis.title=element_blank(),
              axis.text=element_blank(),
              axis.ticks=element_blank(),
              panel.grid=element_blank(),
              plot.margin = unit(c(7,7,7,-5), "pt"))


    gg_ctrls_gtable <- gg_ctrls %>%
        ggplot_build() %>% ggplot_gtable()


    gg_ctrls_gtable$heights <- gg_scatter_gtable$heights

    gg_combined <-
        plot_grid(gg_scatter_gtable,
                  gg_ctrls_gtable,
                  ncol=2, align="none", rel_widths = c(5,1)) +
        draw_plot(gg_rsquareds, x=0.6, y=0.7, width=0.2, height=0.25)

    return(gg_combined)
}
