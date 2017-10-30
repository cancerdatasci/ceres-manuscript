library(ProjectTemplate)
load.project(override.config=list(cache_loading=F, munging=F))


out_dir <- file.path("./output/figures/guide_activity", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


ceres_input_avana <- readRDS("./data/ceres_cache_final/ceres_input_avana.rds")
ceres_output_avana <- readRDS("./data/ceres_cache_final/ceres_output_avana.rds")

ceres_input_gecko <- readRDS("./data/ceres_cache_final/ceres_input_gecko.rds")
ceres_output_gecko <- readRDS("./data/ceres_cache_final/ceres_output_gecko.rds")

ceres_input_wang <- readRDS("./data/ceres_cache_final/ceres_input_wang.rds")
ceres_output_wang <- readRDS("./data/ceres_cache_final/ceres_output_wang.rds")



avana_rs2_df <- read_tsv("./data/raw/avana_rs2.txt") %>%
    select(Guide = `sgRNA Sequence`,
           RS2 = `Ruleset 2 Score`) %>%
    group_by(Guide) %>%
    summarise(RS2 = max(RS2))

gecko_rs2_df <- read_tsv("./data/raw/gecko_rs2.txt") %>%
    select(Guide = `sgRNA Sequence`,
           RS2 = `Ruleset 2 Score`) %>%
    group_by(Guide) %>%
    summarise(RS2 = max(RS2))




### Guide Activity Discrete


activity_breaks <- c(0, 0.2, 0.9, 1)
activity_labels <- c("Low", "Moderate", "High") %>%
    str_c("\n(", activity_breaks %>% head(-1), "-",
          activity_breaks %>% tail(-1), ")")
# activity_colors <- brewer_pal(palette="Set1")(3) %>%
#     set_names(activity_labels)
activity_colors <- viridis::viridis_pal(begin=0, end=0.9, option="D")(3) %>%
    set_names(activity_labels)



library_df <- bind_rows(
    ceres_output_avana$guide_activity %>% mutate(Library="Avana"),
    ceres_output_gecko$guide_activity %>% mutate(Library="GeCKOv2"),
    ceres_output_wang$guide_activity %>% mutate(Library="Wang2017")) %>%
    filter(!is.na(Activity)) %>%
    mutate(ActivityBin = cut(Activity, breaks=activity_breaks,
                          labels=activity_labels),
           Library = factor(Library,
                            levels=c("GeCKOv2", "Avana", "Wang2017"), ordered=T))

library_count_df <- library_df %>% count(Library) %>% mutate(n=prettyNum(n, big.mark=","))

library_df %>% filter(Library %in% c("GeCKOv2", "Avana")) %>%
    ggplot() +
    geom_bar(aes(x=Library, fill=ActivityBin), position="fill", alpha=0.75) +
    geom_text(aes(x=Library, y=1, label=n), vjust=-0.2, size=3,
              data=library_count_df %>% filter(Library %in% c("GeCKOv2", "Avana"))) +
    scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.5, 1), expand=c(0,0)) +
    scale_fill_manual(values=activity_colors) +
    labs(x=NULL, y="Fraction of sgRNAs", title="Guide Activity Score") +
    theme(axis.text.x=element_text(size=12, face="bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.height=unit(2, "lines"),
          legend.margin=margin(0,0,0,0, "pt")) +
    ggsave(file.path(out_dir, "guide_activity_discrete_library.pdf"),
           width=3.5, height=4)

library_df %>% filter(Library == "Wang2017") %>%
    ggplot() +
    geom_bar(aes(x=Library, fill=ActivityBin), position="fill", alpha=0.75) +
    geom_text(aes(x=Library, y=1, label=n), vjust=-0.2, size=3,
              data=library_count_df %>% filter(Library=="Wang2017")) +
    scale_y_continuous(limits=c(0,1.05), breaks=c(0, 0.5, 1), expand=c(0,0)) +
    scale_fill_manual(values=activity_colors) +
    labs(x=NULL, y="Fraction of sgRNAs", title="Guide Activity Score") +
    theme(axis.text.x=element_text(size=12, face="bold"),
          legend.title=element_blank(),
          legend.text=element_text(size=10),
          legend.key.height=unit(2, "lines"),
          legend.margin=margin(0,0,0,0, "pt")) +
    ggsave(file.path(out_dir, "guide_activity_discrete_wang.pdf"),
           width=2.5, height=4)




rs2_breaks <- c(0,0.2,0.4,0.6,0.8,1)
rs2_labels <- c("0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0")
rs2_quantiles <- 5

rs2_df <- library_df %>% inner_join(bind_rows(avana_rs2_df, gecko_rs2_df)) %>%
    filter(RS2 >=0, RS2 <= 1) %>%
    mutate(RS2_Bin = cut(RS2, breaks=rs2_breaks, labels=rs2_labels)) %>%
    group_by(Library) %>%
    mutate(RS2_Quantile = cut(RS2, quantile(RS2, seq(0, 1, by=1/rs2_quantiles)),
                              labels=1:rs2_quantiles, include.lowest=T)) %>%
    ungroup

rs2_count_df <- rs2_df %>%
    count(Library, RS2_Bin) %>%
    ungroup() %>%
    mutate(n=prettyNum(n, big.mark=",")) %>%
    mutate(start = str_extract(RS2_Bin, "^([\\.0-9]+)") %>% as.numeric(),
           end = str_extract(RS2_Bin, "([\\.0-9]+)$") %>% as.numeric(),
           x = (start+end) / 2)
rs2_df %>%
    filter(Library %in% c("Avana", "GeCKOv2")) %>%
    ggplot() +
    facet_wrap(~ Library, scales="free", ncol=2) +
    geom_bar(aes(x=RS2_Quantile, fill=ActivityBin), position="fill", alpha=0.75, show.legend=F) +
    scale_y_continuous(limits=c(0,1), breaks=c(0, 0.5, 1), expand=c(0,0)) +
    labs(x="Doench-Root Score Quintile", y="Fraction of sgRNAs") +
    scale_fill_manual(values=activity_colors) +
    theme(strip.background=element_blank(),
          strip.text=element_text(face="bold", size=12),
          axis.text.x=element_text(size=11, angle=0, vjust=1, hjust=0.5)) +
    ggsave(file.path(out_dir, "guide_activity_discrete_d-r-quartile.pdf"),
           width=7, height=4)

rs2_df %>%
    filter(Library %in% c("Avana", "GeCKOv2")) %>%
    ggplot() +
    facet_wrap(~ Library, scales="free", ncol=2) +
    geom_histogram(aes(x=RS2, fill=ActivityBin),
                   position="fill", alpha=0.75, show.legend=F, breaks=rs2_breaks) +
    geom_vline(xintercept=rs2_breaks, color="white", lwd=2) +
    geom_text(aes(x=x, y=1, label=n), vjust=-0.2, size=3,
              data=rs2_count_df %>% filter(Library %in% c("Avana", "GeCKOv2"))) +
    scale_x_continuous(breaks=rs2_breaks, expand=c(0, 0.025)) +
    scale_y_continuous(limits=c(0, 1.05), breaks=c(0, 0.5, 1), expand=c(0, 0)) +
    scale_fill_manual(values=activity_colors) +
    labs(x="Doench-Root Score", y="Fraction of sgRNAs") +
    theme(strip.background=element_blank(),
          strip.text=element_text(face="bold", size=12),
          panel.spacing=unit(2, "lines")) +
    # plot.margin=unit(c(7,12,7,7), "pt")) +
    # axis.text.x=element_text(size=11, angle=0, vjust=1, hjust=0.5)) +
    ggsave(file.path(out_dir, "guide_activity_discrete_d-r-score.pdf"),
           width=7, height=4)




### Gecko-Avana overlapping guides


# activity_breaks <- c(0, 0.2, 0.9, 1)
# activity_labels <- c("Low", "Moderate", "High") %>%
#     str_c("\n(", activity_breaks %>% head(-1), "-",
#           activity_breaks %>% tail(-1), ")")
# activity_colors <- brewer_pal(palette="Set1")(3) %>%
#     set_names(activity_labels)

avana_gecko_overlap <-
    guide_activity_overlap(ceres_output_gecko$guide_activity,
                           ceres_output_avana$guide_activity,
                           "GeCKOv2", "Avana",
                           activity_breaks=activity_breaks,
                           activity_labels=activity_labels,
                           activity_colors=activity_colors)

ggsave(filename=file.path(out_dir, "guide_overlap_avana_gecko.pdf"),
       avana_gecko_overlap$plot, width=4, height=4)


pdf(file.path(out_dir, "guide_overlap_legend.pdf"), width=1, height=0.5)
grid.draw(avana_gecko_overlap$legend)
dev.off()

avana_wang_overlap <-
    guide_activity_overlap(ceres_output_wang$guide_activity,
                           ceres_output_avana$guide_activity,
                           "Wang2017", "Avana",
                           axis_ticks=c(1, 5000, 10000, 15000),
                           activity_breaks=activity_breaks,
                           activity_labels=activity_labels,
                           activity_colors=activity_colors)

ggsave(filename=file.path(out_dir, "guide_overlap_avana_wang.pdf"),
       avana_wang_overlap$plot, width=4, height=4)


