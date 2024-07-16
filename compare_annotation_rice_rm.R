#!/usr/bin/env Rscript
library(rtracklayer)
library(reshape2)
library(stringr)
library(fmsb)
library(parallel)
library(igraph)
library(ggvenn)
library(ChIPpeakAnno)
library(limma)
# Load custom utility functions for annotation comparisons
source("compare_annot_utils.R")
# standard annot for evaluation:
annot_std_gr <- import("rice_v7/riceTElib/rice_v7.fasta.gff3")
output_dir <- "rice_v7_plots_rm"
dir.create(output_dir, showWarnings = FALSE)

# annotation obtain with RM based on clustom libraries create from Inpactor2, edta and dante_ltr
annot_test_str <- c(
  DANTE_LTR="reference_genomes/rice_v7/libraries/dante_ltr_repre/rice_v7.fasta.gff3",
  Inpactor2 = "reference_genomes/rice_v7/libraries/inpactor2/rice_v7.fasta.gff3",
  EDTA = "reference_genomes/rice_v7/libraries/edta/rice_v7.fasta.gff3"
)

# genome size from ref seq
library(Biostrings)
# Load the genome sequence and calculate sequence lengths
s <- readDNAStringSet("rice_v7/rice_v7.fasta")
SL <- seqlengths(s)
# Read conversion table mapping rice to rexdb names for annotation consistency
rice2rexdb <- read.table("reference_genomes/rice_v7/riceTElib//rice7.0.0.liban_unique_categories_to_rexdb.csv",
  sep = "\t", header = FALSE, as.is = TRUE)
rownames(rice2rexdb) <- rice2rexdb[,1]

# Calculate the total genome size for further calculations
gs <- sum(SL)
# Calculate the total genome size for further calculations
genome_gr <- GRanges(seqnames = names(SL), ranges = IRanges(start = 1, end = SL))
# Calculate the total genome size for further calculations
seqlengths(annot_std_gr) <- SL[seqlevels(annot_std_gr)]
seqlengths(genome_gr) <- SL[seqlevels(genome_gr)]
# Invert intervals in standard annotations and label unannotated regions
annot_str_gr_inv <- setdiff(genome_gr, annot_std_gr, ignore.strand = TRUE)
annot_str_gr_inv$Name <- "no_annotation"
annot_std <-  c(annot_std_gr, annot_str_gr_inv)
annot_std$ori_name <- annot_std$Name
annot_std$Name <- rice2rexdb[annot_std$Name,2]
annot_std$Name[is.na(annot_std$Name)] <- "Uknown"
# Import test annotations using defined paths
annot_test_str_gr_raw <- sapply(annot_test_str, import)
annot_test_str_gr_raw_reduce <- sapply(annot_test_str_gr_raw, reduce)
# calculate overlaps for venn diagram - per base

# we need to compare only True positive
# 1. split groups - LTR, LTR/copia, LTR/gypsy for each method
# 2. split refence annotation also to LTR, LTR/copia, LTR/gypsy
# 3. compare each group separately, for each method keep only overlapping regions

dante_ltr_LTR <- annot_test_str_gr_raw$DANTE_LTR[grepl("LTR$", annot_test_str_gr_raw$DANTE_LTR$Name)]
dante_ltr_Copia <- annot_test_str_gr_raw$DANTE_LTR[grepl("copia", annot_test_str_gr_raw$DANTE_LTR$Name)]
dante_ltr_Gypsy <- annot_test_str_gr_raw$DANTE_LTR[grepl("gypsy", annot_test_str_gr_raw$DANTE_LTR$Name)]

edta_ltr_LTR <- annot_test_str_gr_raw$EDTA[grepl("LTR$", annot_test_str_gr_raw$EDTA$Name) | grepl("LTR/unknown$", annot_test_str_gr_raw$EDTA$Name)]
edta_ltr_Copia <- annot_test_str_gr_raw$EDTA[grepl("Copia", annot_test_str_gr_raw$EDTA$Name)]
edta_ltr_Gypsy <- annot_test_str_gr_raw$EDTA[grepl("Gypsy", annot_test_str_gr_raw$EDTA$Name)]

inpactor_ltr_LTR <- annot_test_str_gr_raw$Inpactor2[grepl("Unknown$", annot_test_str_gr_raw$Inpactor2$Name)]
inpactor_ltr_Copia <- annot_test_str_gr_raw$Inpactor2[grepl("RLC", annot_test_str_gr_raw$Inpactor2$Name)]
inpactor_ltr_Gypsy <- annot_test_str_gr_raw$Inpactor2[grepl("RLG", annot_test_str_gr_raw$Inpactor2$Name)]

annot_std_LTR <- annot_std[grepl("LTR$", annot_std$Name)]
annot_std_Copia <- annot_std[grepl("copia", annot_std$Name)]
annot_std_Gypsy <- annot_std[grepl("gypsy", annot_std$Name)]

# get True Positive - this is like bedtools intersect
dante_ltr_LTR_TP <- intersect(dante_ltr_LTR, annot_std_LTR)
dante_ltr_Copia_TP <- intersect(dante_ltr_Copia, annot_std_Copia)
dante_ltr_Gypsy_TP <- intersect(dante_ltr_Gypsy, annot_std_Gypsy)

edta_ltr_LTR_TP <- intersect(edta_ltr_LTR, annot_std_LTR)
edta_ltr_Copia_TP <- intersect(edta_ltr_Copia, annot_std_Copia)
edta_ltr_Gypsy_TP <- intersect(edta_ltr_Gypsy, annot_std_Gypsy)

inpactor_ltr_LTR_TP <- intersect(inpactor_ltr_LTR, annot_std_LTR)
inpactor_ltr_Copia_TP <- intersect(inpactor_ltr_Copia, annot_std_Copia)
inpactor_ltr_Gypsy_TP <- intersect(inpactor_ltr_Gypsy, annot_std_Gypsy)

# get False Positive - by subtracting TP from all

dante_ltr_FP <- subtract(
  c(annot_test_str_gr_raw_reduce$DANTE_LTR),
  c(dante_ltr_LTR_TP, dante_ltr_Copia_TP, dante_ltr_Gypsy_TP)
) |> unlist()

inpactor_ltr_FP <- subtract(
  c(annot_test_str_gr_raw_reduce$Inpactor2),
  c(inpactor_ltr_LTR_TP, inpactor_ltr_Copia_TP, inpactor_ltr_Gypsy_TP)
) |> unlist()

edta_ltr_FP <- subtract(
  c(annot_test_str_gr_raw_reduce$EDTA),
  c(edta_ltr_LTR_TP, edta_ltr_Copia_TP, edta_ltr_Gypsy_TP)
) |> unlist()



annot_test_str_gr_raw_TP <- list(DANTE_LTR=c(dante_ltr_LTR_TP, dante_ltr_Copia_TP, dante_ltr_Gypsy_TP),
                                 Inpactor2=c(inpactor_ltr_LTR_TP, inpactor_ltr_Copia_TP, inpactor_ltr_Gypsy_TP),
                                 EDTA=c(edta_ltr_LTR_TP, edta_ltr_Copia_TP, edta_ltr_Gypsy_TP))

annot_test_str_gr_raw_FP <- lapply(list(DANTE_LTR=dante_ltr_FP,
                                 Inpactor2=inpactor_ltr_FP,
                                 EDTA=edta_ltr_FP), reduce)


annot_test_str_gr_raw_TP_reduce <- sapply(annot_test_str_gr_raw_TP, reduce)
res <- makeVennDiagram(annot_test_str_gr_raw_TP_reduce, by='base', NameOfPeaks = names(annot_test_str_gr_raw_TP_reduce))
VennCountsDF_TP <- as.data.frame.matrix(res$vennCounts)
VennCountsDF_TP$Perc <- VennCountsDF_TP$Count / sum(VennCountsDF_TP$Count) * 100
write.table(VennCountsDF_TP, paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM_TP.csv"), sep = "\t", quote = FALSE, row.names = FALSE)


res_FP <- makeVennDiagram(annot_test_str_gr_raw_FP, by='base', NameOfPeaks = names(annot_test_str_gr_raw_FP))
VennCountsDF_FP <- as.data.frame.matrix(res_FP$vennCounts)
VennCountsDF_FP$Perc <- VennCountsDF_FP$Count / sum(VennCountsDF_FP$Count) * 100
write.table(VennCountsDF_FP, paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM_FP.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

annot_test_str_gr_raw_TP_reduce <- sapply(annot_test_str_gr_raw_TP, reduce)
annot_test_str_gr_raw_FP_reduce <- sapply(annot_test_str_gr_raw_FP, reduce)


FP_size <- sapply(annot_test_str_gr_raw_FP_reduce, function(x)sum(width(x)))
TP_size <- sapply(annot_test_str_gr_raw_TP_reduce, function(x)sum(width(x)))
# stack barplot TP on top of FP

pdf(paste0(output_dir, "/annot_stat_comparison_maize_TP_FP_RM.pdf"), width = 3, height = 5, pointsize = 3)
par(mar = c(10, 7, 1.5, 1), lwd=3, mgp = c(3.5, 0.2, 0))
barplot(rbind(TP_size/1e6, FP_size/1e6), beside = FALSE, col = c("#AAFFAA", "#FFAAAA"),
        names.arg = names(FP_size), las = 2, cex.names = 2, cex.axis = 1.5, cex.lab = 2,
        ylab = "Size (Mbp)", xpd = FALSE, axes = FALSE, ylim = c(0,120), xlim =c(0,6))
par(mgp = c(3, 1, 0))
axis(2, lwd = 2, cex.lab = 4, cex.axis = 2, at = seq(0, 100, 20), labels = seq(0, 100, 20))
#legend("top", legend = c("True positive", "False positive"), fill = c("#AAFFAA", "#FFAAAA"), cex = 2.5, bty = "n")
dev.off()
save.image(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))





pdf(paste0(output_dir, "/annotation_overlaps_per_base_venn_RM.pdf"), width = 13, height = 13)
res <- makeVennDiagram(annot_test_str_gr_raw_reduce, by='base', NameOfPeaks = names(annot_test_str_gr_raw_reduce))
dev.off()
# export counts
VennCountsDF <- as.data.frame.matrix(res$vennCounts)
VennCountsDF$Perc <- VennCountsDF$Count / sum(VennCountsDF$Count) * 100
write.table(VennCountsDF, paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.csv"), sep = "\t", quote = FALSE, row.names = FALSE)


for (i in names(annot_test_str_gr_raw)){
  seqlengths(annot_test_str_gr_raw[[i]]) <- SL[seqlevels(annot_test_str_gr_raw[[i]])]
}


annot_intact_complete <- sapply(annot_test_str_gr_raw, make_complete, genome_gr)
annot_intact_complete$DANTE_LTR$source <- 'dante_ltr'
annot_intact_complete$Inpactor2$source <- 'Inpactor2'
annot_intact_complete$EDTA$source <- 'EDTA'
annot_intact_complete_adj <- sapply(annot_intact_complete, adjust_names2)

annot_pairs <- lapply(annot_intact_complete_adj, extract_all_matches,annot_std)
annot_confusion_matrix <- lapply(annot_pairs, get_confusion_matrix)
annot_stat <- lapply(annot_confusion_matrix, get_stat_from_confusion_matrix)

annot_pairs_gr_list <- lapply(annot_intact_complete_adj, get_matched_granges,annot_std)

# comparison of all LTR-RT elements, for the purpose of comparison, copia annotated as gypsy is considered false positive
stat_LTR_RT <- lapply(annot_pairs, calculate_TP_TF_FP_FN_from_pairs_LTR_category)
stat_values_LTR_RT <- lapply(stat_LTR_RT, calculate_statistics_from_groups)
stat_values_LTR_RT <- as.data.frame.matrix(do.call(rbind, stat_values_LTR_RT))
stat_values_LTR_RT_rc <-  rbind(rep(1, 6), rep(0,6), stat_values_LTR_RT)
# this is to set the scale

stat_to_export <- cbind(do.call(rbind, stat_LTR_RT), stat_values_LTR_RT)
stat_to_export$Method <- names(stat_LTR_RT)
write.table(stat_to_export, paste0(output_dir, "/stat_values_LTR_RT.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(paste0(output_dir, "/annot_stat_comparison_full_elements_radarchart_ltr_RM.pdf"), width = 10, height = 5, pointsize = 3)
par(mfrow = c(1,2))
vennDiagram(res$vennCounts, lwd=2, cex = 1)
colors <- c(DANTE_LTR = "red", Inpactor2 = "blue", EDTA = "green")
radarchart(stat_values_LTR_RT_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=3,
           title="", plwd = 4, pcol = c("red", "blue", "green", "black", "orange", "purple"))
# do not draw
legend(-1.4,-0.8, legend = names(colors), col = colors, lwd = 3, cex = 1.5, bty = 'n')
dev.off()




# export annot pairs as csv table
write.table(annot_pairs$DANTE_LTR, paste0(output_dir, "/annot_pairs_DANTE_LTR.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(annot_pairs$Inpactor2, paste0(output_dir, "/annot_pairs_Inpactor2.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(annot_pairs$EDTA, paste0(output_dir, "/annot_pairs_EDTA.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




annot_intact_complete_ajd_ltr_only <- lapply(annot_intact_complete_adj, function(x) {
  x$Name <- gsub("LTR[|].*", "LTR",  x$Name)
  x$Name <- gsub("^LTR$", "Class_I|LTR", x$Name)
  x
})



annot_str_ltr_only <- annot_std
annot_str_ltr_only$Name <- gsub("LTR[|].*", "LTR",  annot_std$Name)

annot_pairs_ltr_only <- lapply(annot_intact_complete_ajd_ltr_only, extract_all_matches,annot_str_ltr_only)
annot_confusion_matrix_ltr_only<- lapply(annot_pairs_ltr_only, get_confusion_matrix)
annot_stat_ltr_only <- lapply(annot_confusion_matrix_ltr_only, get_stat_from_confusion_matrix)


# export tables from annot_stat to csv
annot_stat_df <- as.data.frame(do.call(rbind, annot_stat))
# include only copia and gypsy rows and reorder columns
cols <- c('sensitivity', 'specificity', 'accuracy', 'F1', 'precision', 'FDR',
          'TP', 'TN', 'FP', 'FN')
annot_stat_df <- annot_stat_df[grepl("copia|gypsy", rownames(annot_stat_df)),cols]

df2 <- as.data.frame.matrix(do.call(rbind, strsplit(rownames(annot_stat_df), "[.]")))
colnames(df2) <- c("method", "class")

annot_stat_df2 <- cbind(df2, annot_stat_df)

write.table(annot_stat_df2, file = paste0(output_dir, "/annot_stat_using_RM_annotated_genome.csv")
  , sep = "\t", quote = FALSE, row.names = FALSE)


# export annotation pairs

write.table(annot_pairs$DANTE_LTR, file = paste0(output_dir, "/annot_pairs_using_RM_annotated_genome_DANTE_LTR.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annot_pairs$EDTA, file = paste0(output_dir,"/annot_pairs_using_RM_annotated_genome_EDTA.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(annot_pairs$Inpactor2, file = paste0(output_dir, "/annot_pairs_using_RM_annotated_genome_Inpactor2.csv"), sep = "\t", quote = FALSE, row.names = FALSE)


save.image(paste0(output_dir, "/annot_stat_using_RM_annotated_genome.RData"))
