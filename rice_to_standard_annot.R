#!/usr/bin/env Rscript
library(rtracklayer)
source("reference_genomes/compare_annot_utils.R")
annot_test_str <- c(
  DANTE_LTR="reference_genomes/rice_v7/rice_v7.fasta_dante_ltr.gff3",
  Inpactor2 = "reference_genomes/rice_v7/rice_v7.fasta_inpactor2_c3/Inpactor2_predictions.bed",
  EDTA = "reference_genomes/rice_v7/rice_v7.fasta.mod.EDTA.raw/rice_v7.fasta.mod.LTR.intact.gff3"
)
annot_test_str_gr_raw <- sapply(annot_test_str, import)

dante_ltr <- annot_test_str_gr_raw$DANTE_LTR[annot_test_str_gr_raw$DANTE_LTR$type == "transposable_element",]
dante_ltr <- dante_ltr[dante_ltr$Rank != "D"]
dante_ltr$Name[grepl('gypsy', dante_ltr$Name)] <- "LTR/Gypsy"
dante_ltr$Name[grepl('copia', dante_ltr$Name)] <- "LTR/Copia"
dante_ltr <- remove_overlaping_elements(dante_ltr)

inpactor2 <- annot_test_str_gr_raw$Inpactor2
inpactor2$Name <- inpactor2$name
inpactor2$Name[grepl('RLG', inpactor2$Name)] <- "LTR/Gypsy"
inpactor2$Name[grepl('RLC', inpactor2$Name)] <- "LTR/Copia"
inpactor2 <- remove_overlaping_elements(inpactor2)

edta <- annot_test_str_gr_raw$EDTA[annot_test_str_gr_raw$EDTA$type == "repeat_region"]
edta$Name <- edta$Classification
edta <- remove_overlaping_elements(edta)

# export to gff3
export(dante_ltr, "reference_genomes/rice_v7/rice_v7.fasta_dante_ltr.standardised.gff3", format = "gff3")
export(inpactor2, "reference_genomes/rice_v7/rice_v7.fasta_inpactor2_c3/Inpactor2_predictions_standartised.gff3", format = "gff3")
export(edta,  "reference_genomes/rice_v7/rice_v7.fasta.mod.EDTA.raw/rice_v7.fasta.mod.LTR.intact_standartised.gff3", format = "gff3")
