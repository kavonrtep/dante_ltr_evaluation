library(rtracklayer)
source("compare_annot_utils.R")
# adjust reference annotation
annot_std_gr <- import("reference_genomes/maize_B73/v3_MIPS_Repeats_loci_adjusted_names.gff3")

annot_std_gr_standardised <- annot_std_gr
annot_std_gr_standardised$Name[annot_std_gr_standardised$Name == "Class_I|LTR|Ty1/copia"] <- "LTR/Copia"
annot_std_gr_standardised$Name[annot_std_gr_standardised$Name == "Class_I|LTR|Ty3/gypsy"] <- "LTR/Gypsy"
annot_std_gr_standardised$Name[annot_std_gr_standardised$Name == "Class_I|LTR|"] <- "LTR"

annt_test_str <- c(
  DANTE_LTR="reference_genomes/maize_B73/B73_RefGen_v3.fa_dante_ltr.gff3",
  Inpactor2 = "reference_genomes/maize_B73/B73_RefGen_v3.fa_inpactor2_c3/Inpactor2_predictions.bed",
  EDTA = "reference_genomes/maize_B73/B73_RefGen_v3.fa.mod.EDTA.raw/B73_RefGen_v3.fa.mod.LTR.intact.gff3"
)
annot_test_str_gr_raw <- sapply(annt_test_str, import)

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
export(dante_ltr, "reference_genomes/maize_B73/B73_RefGen_v3.fa_dante_ltr.standardised.gff3", format = "gff3")
export(inpactor2, "reference_genomes/maize_B73/B73_RefGen_v3.fa_inpactor2_c3/Inpactor2_predictions_standartised.gff3", format = "gff3")
export(edta,  "reference_genomes/maize_B73/B73_RefGen_v3.fa.mod.EDTA.raw/B73_RefGen_v3.fa.mod.LTR.intact_standartised.gff3", format = "gff3")

# export standard annot
export(annot_std_gr_standardised, "reference_genomes/maize_B73/v3_MIPS_Repeats_loci_adjusted_names_standardised.gff3", format = "gff3")

library(parallel)
# it seems alot of feature is overlaping - this must be handed

annot_std_gr_standardised_clean <- gff_cleanup(annot_std_gr_standardised)

export(annot_std_gr_standardised_clean, "reference_genomes/maize_B73/v3_MIPS_Repeats_loci_adjusted_names_standardised_clean.gff3", format = "gff3")

