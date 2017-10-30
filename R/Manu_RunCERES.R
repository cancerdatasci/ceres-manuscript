library(ProjectTemplate)
load.project()

cn_seg_file <- "./data/raw/CCLE_copynumber_2013-12-03.seg.txt"
gene_annot_file <- "./data/raw/CCDS.current.txt"


### Run CERES on Wang2017 data

wang_dep_file <- "./data/interim/wang_guide_dep.gct"
write.gct(wang_guide_dep, wang_dep_file)

wang_rep_map_file <- "./data/raw/wang_replicate_map.tsv"

wang_inputs_dir <- file.path("./data/ceres_inputs/wang")
dir.create(wang_inputs_dir, recursive=T, showWarnings=F)
wang_outputs_dir <- file.path("./data/ceres_outputs/wang")
dir.create(wang_outputs_dir, recursive=T, showWarnings=F)

prepare_ceres_inputs(inputs_dir=wang_inputs_dir,
                     dep_file=wang_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=wang_rep_map_file,
                     genome_id="hg19",
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

wang_ceres <-
    wrap_ceres(sg_path=file.path(wang_inputs_dir, "guide_sample_dep.rds"),
               cn_path=file.path(wang_inputs_dir, "locus_sample_cn.rds"),
               guide_locus_path=file.path(wang_inputs_dir, "guide_locus.rds"),
               locus_gene_path=file.path(wang_inputs_dir, "locus_gene.rds"),
               replicate_map_path=file.path(wang_inputs_dir, "replicate_map.rds"),
               run_id="wang2017",
               params=list(lambda_g=0.68129207))

wang_ceres_scaled <-
    scale_to_essentials(wang_ceres$gene_essentiality_results$ge_fit)

saveRDS(wang_ceres_scaled,
        file.path(wang_outputs_dir, "ge_fit.Rds"))
saveRDS(wang_ceres$cutting_effect_results$ce_df,
        file.path(wang_outputs_dir, "ce_fit.Rds"))
saveRDS(wang_ceres$sgrna_results$sgrna_df,
        file.path(wang_outputs_dir, "sgRNA_fit.Rds"))


### Run CERES on GeCKOv2 data

gecko_dep_file <- "./data/interim/gecko_guide_dep.gct"
write.gct(gecko_guide_dep, gecko_dep_file)

gecko_rep_map_file <- "./data/raw/gecko_replicate_map.tsv"

gecko_inputs_dir <- file.path("./data/ceres_inputs/gecko")
dir.create(gecko_inputs_dir, recursive=T, showWarnings=F)
gecko_outputs_dir <- file.path("./data/ceres_outputs/gecko")
dir.create(gecko_outputs_dir, recursive=T, showWarnings=F)

prepare_ceres_inputs(inputs_dir=gecko_inputs_dir,
                     dep_file=gecko_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=gecko_rep_map_file,
                     genome_id="hg19",
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

gecko_ceres <-
    wrap_ceres(sg_path=file.path(gecko_inputs_dir, "guide_sample_dep.rds"),
               cn_path=file.path(gecko_inputs_dir, "locus_sample_cn.rds"),
               guide_locus_path=file.path(gecko_inputs_dir, "guide_locus.rds"),
               locus_gene_path=file.path(gecko_inputs_dir, "locus_gene.rds"),
               replicate_map_path=file.path(gecko_inputs_dir, "replicate_map.rds"),
               run_id="gecko",
               params=list(lambda_g=0.68129207))

gecko_ceres_scaled <-
    scale_to_essentials(gecko_ceres$gene_essentiality_results$ge_fit)

saveRDS(gecko_ceres_scaled,
        file.path(gecko_outputs_dir, "ge_fit.Rds"))
saveRDS(gecko_ceres$cutting_effect_results$ce_df,
        file.path(gecko_outputs_dir, "ce_fit.Rds"))
saveRDS(gecko_ceres$sgrna_results$sgrna_df,
        file.path(gecko_outputs_dir, "sgRNA_fit.Rds"))


### Run CERES on Avana data

avana_dep_file <- "./data/interim/avana_guide_dep.gct"
write.gct(avana_guide_dep, avana_dep_file)

avana_rep_map_file <- "./data/raw/avana_replicate_map.tsv"

avana_inputs_dir <- file.path("./data/ceres_inputs/avana")
dir.create(avana_inputs_dir, recursive=T, showWarnings=F)
avana_outputs_dir <- file.path("./data/ceres_outputs/avana")
dir.create(avana_outputs_dir, recursive=T, showWarnings=F)

prepare_ceres_inputs(inputs_dir=avana_inputs_dir,
                     dep_file=avana_dep_file,
                     cn_seg_file=cn_seg_file,
                     gene_annot_file=gene_annot_file,
                     rep_map_file=avana_rep_map_file,
                     genome_id="hg19",
                     chromosomes=paste0("chr", 1:22),
                     dep_normalize="zmad")

avana_ceres <-
    wrap_ceres(sg_path=file.path(avana_inputs_dir, "guide_sample_dep.rds"),
               cn_path=file.path(avana_inputs_dir, "locus_sample_cn.rds"),
               guide_locus_path=file.path(avana_inputs_dir, "guide_locus.rds"),
               locus_gene_path=file.path(avana_inputs_dir, "locus_gene.rds"),
               replicate_map_path=file.path(avana_inputs_dir, "replicate_map.rds"),
               run_id="avana",
               params=list(lambda_g=0.56234133))

avana_ceres_scaled <-
    scale_to_essentials(avana_ceres$gene_essentiality_results$ge_fit)

saveRDS(avana_ceres_scaled,
        file.path(avana_outputs_dir, "ge_fit.Rds"))
saveRDS(avana_ceres$cutting_effect_results$ce_df,
        file.path(avana_outputs_dir, "ce_fit.Rds"))
saveRDS(avana_ceres$sgrna_results$sgrna_df,
        file.path(avana_outputs_dir, "sgRNA_fit.Rds"))
