
guide_activity_overlap <- function(activity_x, activity_y,
                               label_x="X", label_y="Y",
                               activity_breaks=c(0, 0.2, 0.9, 1),
                               activity_labels=c("Low", "Moderate", "High"),
                               activity_colors=brewer_pal(palette="Set1")(3),
                               axis_ticks=c(1, seq(1000, 5000, by=1000))) {



    label_x <- str_c(label_x, " Guide Activty Rank")
    label_y <- str_c(label_y, " Guide Activty Rank")


    activity_colors <- activity_colors %>%
        set_names(activity_labels)

    activity_df <- inner_join(activity_x, activity_y, by="Guide") %>%
        distinct(Guide, Activity.x, Activity.y) %>%
        filter(!is.na(Activity.x), !is.na(Activity.y)) %>%
        mutate(Rank.x = rank(-round(Activity.x, digits=10), ties.method="random"),
               Rank.y = rank(-round(Activity.y, digits=10), ties.method="random")) %>%
        mutate(Bin.x = cut(Activity.x, activity_breaks, activity_labels),
               Bin.y = cut(Activity.y, activity_breaks, activity_labels))

    cor_results <- cor.test(activity_df$Activity.x, activity_df$Activity.y, method="spearman")


    rank_density <- MASS::kde2d(activity_df$Rank.x, activity_df$Rank.y, n=50)
    rank_density_df <- rank_density$z %>%
        magrittr::set_colnames(rank_density$y) %>%
        magrittr::set_rownames(rank_density$x) %>%
        reshape2::melt() %>%
        magrittr::set_colnames(c("Rank.x", "Rank.y", "Density"))




    rank_dens_gg <-
        ggplot(rank_density_df, aes(Rank.x, Rank.y, fill = Density)) +
        geom_raster() +
        scale_fill_gradientn(colours=scales::brewer_pal(palette="Blues")(8),
                             name="Density",
                             guide=guide_colorbar(direction="horizontal",
                                                  title.position="top",
                                                  title.vjust=0)) +
        scale_y_continuous(expand=c(0,0), breaks=c(1, seq(1000, nrow(activity_df), by=1000))) +
        scale_x_continuous(expand=c(0,0), breaks=c(1, seq(1000, nrow(activity_df), by=1000))) +
        labs(title = str_c(prettyNum(nrow(activity_df),big.mark=","), " Shared sgRNAs"),
             x=label_x,
             y=label_y) +
        theme(panel.border=element_rect(size=1, color="black", linetype=1, fill=NA),
              axis.line.x=element_blank(),
              axis.line.y=element_blank(),
              axis.ticks=element_blank(),
              axis.text=element_blank(),
              axis.title=element_blank(),
              legend.text=element_blank(),
              legend.title=element_text(size=10),
              plot.margin=margin(7, 14, 0, 0, "pt"))
    rank_dens_gtable <- {rank_dens_gg  + theme(legend.position="none")} %>%
        ggplot_build %>% ggplot_gtable

    y_shared_gg <-
        ggplot(activity_df, aes(x=Rank.y, fill=Bin.y)) +
        geom_area(aes(y=1), alpha=0.75, show.legend=F) +
        scale_fill_manual(values=activity_colors) +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0),
                           breaks=axis_ticks) +
        labs(x=label_y) +
        coord_flip() +
        theme(axis.line=element_blank(),
              axis.ticks.x=element_blank(),
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_text(angle=90, hjust=0.5),
              plot.margin=margin(7, -3, 0, 7, "pt"))
    y_shared_gtable <- y_shared_gg %>% ggplot_build %>% ggplot_gtable
    rank_dens_gtable$heights <- unit.pmax(rank_dens_gtable$heights, y_shared_gtable$heights)
    y_shared_gtable$heights <- unit.pmax(rank_dens_gtable$heights, y_shared_gtable$heights)

    x_shared_gg <-
        ggplot(activity_df, aes(x=Rank.x, fill=Bin.x)) +
        geom_area(aes(y=1), alpha=0.75, show.legend=F) +
        scale_fill_manual(values=activity_colors) +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(expand=c(0,0),
                           breaks=axis_ticks) +
        labs(x=label_x) +
        theme(axis.line=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              plot.margin=margin(-3, 14, 7, 0, "pt"))
    x_shared_gtable <- x_shared_gg %>% ggplot_build %>% ggplot_gtable()
    rank_dens_gtable$widths <- unit.pmax(rank_dens_gtable$widths, x_shared_gtable$widths)
    x_shared_gtable$widths <- unit.pmax(rank_dens_gtable$widths, x_shared_gtable$widths)

    null_gg <- ggplot(data_frame()) +
        theme(axis.text=element_blank(),
              axis.line=element_blank(),
              axis.title=element_blank(),
              axis.ticks=element_blank())


    grid.newpage()

    overlap_grob <-
        gridExtra::arrangeGrob(y_shared_gtable, rank_dens_gtable,
                               null_gg, x_shared_gtable,
                               nrow=2, ncol=2,
                               widths=c(1, 4), heights=c(4,1))
    grid.draw(overlap_grob)


    g_legend <- function(a.gplot){
        tmp <- ggplot_gtable(ggplot_build(a.gplot))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)}

    my_legend <- g_legend(rank_dens_gg + theme(legend.title=element_blank()))

    return(list(plot=overlap_grob, legend=my_legend, cor=cor_results))
}
