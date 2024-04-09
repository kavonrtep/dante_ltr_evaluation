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
source("compare_annot_utils.R")
# standard annot for evaluation:
output_dir <- "maize_B73_plots_rm"
dir.create(output_dir, showWarnings = FALSE)

## INPUT DATA
annot_std_gr <- import("/mnt/raid/454_data/dante/reference_genomes/maize_B73/v3_MIPS_Repeats_loci_adjusted_names.gff3")
p <- findOverlaps(annot_std_gr, annot_std_gr)
ovlp <- p[!from(p) == to(p)]
# not empty - the are overlaps!!, must be cleaned up
# remove conflicts
annot_std_gr <- gff_cleanup(annot_std_gr) # this takes ~ 10-15 min
check_for_overlaps(annot_std_gr)
# TODO
annot_test_str <- c(
    DANTE_LTR="/mnt/raid/454_data/dante/reference_genomes/maize_B73/libraries2/dante_ltr/B73_RefGen_v3_noAmbiguities.fa.gff3",
    Inpactor2 = "/mnt/raid/454_data/dante/reference_genomes/maize_B73/libraries2/inpactor/B73_RefGen_v3_noAmbiguities.fa.gff3",
    EDTA = "/mnt/raid/454_data/dante/reference_genomes/maize_B73/libraries2/edta/B73_RefGen_v3_noAmbiguities.fa.gff3"
)
library(Biostrings)
s <- readDNAStringSet("maize_B73/B73_RefGen_v3_noAmbiguities.fa")
names(s) <- gsub(" .+$", "", names(s))
SL <- seqlengths(s)



# gneome size:
gs <- sum(SL)
# whole genome as GRanges
genome_gr <- GRanges(seqnames = names(SL), ranges = IRanges(start = 1, end = SL))
# add seq sizes to all annots
seqlengths(annot_std_gr) <- SL[seqlevels(annot_std_gr)]
seqlengths(genome_gr) <- SL[seqlevels(genome_gr)]
# invert intervals in annot_std_gr and adjust names
annot_str_gr_inv <- setdiff(genome_gr, annot_std_gr, ignore.strand = TRUE)
annot_str_gr_inv$Name <- "no_annotation"
annot_std <-  c(annot_std_gr, annot_str_gr_inv)
annot_std$ori_name <- annot_std$Name

annot_test_str_gr_raw <- sapply(annot_test_str, import)
# calculate overlaps for venn diagram - per base
annot_test_str_gr_raw_reduce <- sapply(annot_test_str_gr_raw, reduce)
save.image(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))
# load(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))



makeVennDiagram_large <- function (gr_list){
    # gr must be reducee!
    labels <- names(gr_list)
    # add labels to the source
    for (i in seq_along(gr_list)){
        gr_list[[i]]$source <- factor(labels[i], levels = labels)
    }
    # merge all
    gr_merged <- unlist(GRangesList(gr_list))
     name_for_split <- as.character(seqnames(gr_merged))
    name_for_split[!grepl("Chr", name_for_split)] <- "Other"
    gr_merged_split <- split(gr_merged, name_for_split)
    res_list <- lapply(gr_merged_split, function(x) {
        xx <- split(x, x$source)
        print("-------")
        print(seqnames(x))
        print("-------")
        res <- makeVennDiagram(xx, by='base', NameOfPeaks = names(xx), plot = FALSE)
        return(res)
    })
    res_vennCounts <- Reduce(function(x, y) {
        x$vennCounts <- x$vennCounts + y$vennCounts
        return(x)
    }, res_list)
    return(res_vennCounts)
}

res_vennCounts <- makeVennDiagram_large(annot_test_str_gr_raw_reduce)
# export counts - thi is all - TP and FP mixed
VennCountsDF <- as.data.frame.matrix(res_vennCounts$vennCounts)
VennCountsDF$Perc <- VennCountsDF$Count / sum(VennCountsDF$Count) * 100
write.table(VennCountsDF, paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
save.image(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))



# we need to compare only True positive and then False positive
# 1. split groups - LTR, LTR/copia, LTR/gypsy for each method
# 2. split refence annotation also to LTR, LTR/copia, LTR/gypsy
# 3. compare each group separately, for each method keep only overlapping regions
dante_ltr_LTR <- annot_test_str_gr_raw$DANTE_LTR[grepl("LTR$", annot_test_str_gr_raw$DANTE_LTR$Name)]
dante_ltr_Copia <- annot_test_str_gr_raw$DANTE_LTR[grepl("Copia", annot_test_str_gr_raw$DANTE_LTR$Name)]
dante_ltr_Gypsy <- annot_test_str_gr_raw$DANTE_LTR[grepl("Gypsy", annot_test_str_gr_raw$DANTE_LTR$Name)]

edta_ltr_LTR <- annot_test_str_gr_raw$EDTA[grepl("LTR$", annot_test_str_gr_raw$EDTA$Name) | grepl("LTR/unknown$", annot_test_str_gr_raw$EDTA$Name)]
edta_ltr_Copia <- annot_test_str_gr_raw$EDTA[grepl("Copia", annot_test_str_gr_raw$EDTA$Name)]
edta_ltr_Gypsy <- annot_test_str_gr_raw$EDTA[grepl("Gypsy", annot_test_str_gr_raw$EDTA$Name)]

inpactor_ltr_LTR <- annot_test_str_gr_raw$Inpactor2[grepl("LTR$$", annot_test_str_gr_raw$Inpactor2$Name)]
inpactor_ltr_Copia <- annot_test_str_gr_raw$Inpactor2[grepl("Copia", annot_test_str_gr_raw$Inpactor2$Name)]
inpactor_ltr_Gypsy <- annot_test_str_gr_raw$Inpactor2[grepl("Gypsy", annot_test_str_gr_raw$Inpactor2$Name)]

annot_std_LTR <- annot_std[grepl("LTR$", annot_std$Name)]
annot_std_Copia <- annot_std[grepl("copia", annot_std$Name)]
annot_std_Gypsy <- annot_std[grepl("gypsy", annot_std$Name)]


# get True Positive - this is like bedtools intersect
dante_ltr_LTR_TP <- GenomicRanges::intersect(dante_ltr_LTR, annot_std_LTR)
dante_ltr_Copia_TP <- GenomicRanges::intersect(dante_ltr_Copia, annot_std_Copia)
dante_ltr_Gypsy_TP <- GenomicRanges::intersect(dante_ltr_Gypsy, annot_std_Gypsy)

edta_ltr_LTR_TP <- GenomicRanges::intersect(edta_ltr_LTR, annot_std_LTR)
edta_ltr_Copia_TP <- GenomicRanges::intersect(edta_ltr_Copia, annot_std_Copia)
edta_ltr_Gypsy_TP <- GenomicRanges::intersect(edta_ltr_Gypsy, annot_std_Gypsy)

inpactor_ltr_LTR_TP <- GenomicRanges::intersect(inpactor_ltr_LTR, annot_std_LTR)
inpactor_ltr_Copia_TP <- GenomicRanges::intersect(inpactor_ltr_Copia, annot_std_Copia)
inpactor_ltr_Gypsy_TP <- GenomicRanges::intersect(inpactor_ltr_Gypsy, annot_std_Gypsy)

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
res_tp <- makeVennDiagram_large(annot_test_str_gr_raw_TP_reduce)
VennCountsDF_TP <- as.data.frame.matrix(res_tp$vennCounts)
VennCountsDF_TP$Perc <- VennCountsDF_TP$Count / sum(VennCountsDF_TP$Count) * 100
write.table(VennCountsDF_TP, paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM_TP.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

annot_test_str_gr_raw_FP_reduce <- sapply(annot_test_str_gr_raw_FP, reduce)
res_fp <- makeVennDiagram_large(annot_test_str_gr_raw_FP_reduce)
VennCountsDF_FP <- as.data.frame.matrix(res_fp$vennCounts)
VennCountsDF_FP$Perc <- VennCountsDF_FP$Count / sum(VennCountsDF_FP$Count) * 100
write.table(VennCountsDF_FP, paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM_FP.csv"), sep = "\t", quote = FALSE, row.names = FALSE)



# plot barplot - FP/TR
FP_size <- sapply(annot_test_str_gr_raw_FP_reduce, function(x)sum(width(x)))
TP_size <- sapply(annot_test_str_gr_raw_TP_reduce, function(x)sum(width(x)))
# stack barplot TP on top of FP

pdf(paste0(output_dir, "/annot_stat_comparison_maize_TP_FP_RM.pdf"), width = 3, height = 5, pointsize = 3)
par(mar = c(10, 7, 1.5, 1), lwd=3, mgp = c(3.5, 0.2, 0))
barplot(rbind(TP_size/1e9, FP_size/1e9), beside = FALSE, col = c("#AAFFAA", "#FFAAAA"),
        names.arg = names(FP_size), las = 2, cex.names = 2, cex.axis = 1.5, cex.lab = 2,
        ylab = "Size (Gbp)", xpd = FALSE, axes = FALSE, ylim = c(0,1.99), xlim =c(0,6))
par(mgp = c(3, 1, 0))
axis(2, lwd = 2, cex.lab = 4, cex.axis = 2)
legend("top", legend = c("True positive", "False positive"), fill = c("#AAFFAA", "#FFAAAA"), cex = 2.5, bty = "n")
dev.off()
save.image(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))

# load(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))
annot_intact_complete <- sapply(annot_test_str_gr_raw, make_complete, genome_gr)
annot_intact_complete$DANTE_LTR$source <- 'dante_ltr'
annot_intact_complete$Inpactor2$source <- 'Inpactor2'
annot_intact_complete$EDTA$source <- 'EDTA'
annot_intact_complete_adj <- sapply(annot_intact_complete, adjust_names2)

annot_pairs <- lapply(annot_intact_complete_adj, extract_all_matches,annot_std)
annot_confusion_matrix <- lapply(annot_pairs, get_confusion_matrix)
annot_stat <- lapply(annot_confusion_matrix, get_stat_from_confusion_matrix)

# comparison of all LTR-RT elements, for the purpose of comparison, copia annotated as gypsy is considered false positive
stat_LTR_RT <- lapply(annot_pairs, calculate_TP_TF_FP_FN_from_pairs_LTR_category)
stat_values_LTR_RT <- lapply(stat_LTR_RT, calculate_statistics_from_groups)
stat_values_LTR_RT <- as.data.frame.matrix(do.call(rbind, stat_values_LTR_RT))
stat_values_LTR_RT_rc <-  rbind(rep(1, 6), rep(0,6), stat_values_LTR_RT) # this is to set the scale
stat_to_export <- cbind(do.call(rbind, stat_LTR_RT), stat_values_LTR_RT)
stat_to_export$Method <- names(stat_LTR_RT)



pdf(paste0(output_dir, "/annot_stat_comparison_full_elements_radarchart_ltr_RM.pdf"), width = 10, height = 5, pointsize = 3)
par(mfrow = c(1,2))
vennDiagram(res_vennCounts$vennCounts, lwd=2, cex = 2)
colors <- c(DANTE_LTR = "red", Inpactor2 = "blue", EDTA = "green")
radarchart(stat_values_LTR_RT_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=3,
           title="", plwd = 4, pcol = c("red", "blue", "green", "black", "orange", "purple"))
# do not draw
legend(-1.4,-0.8, legend = names(colors), col = colors, lwd = 3, cex = 1.5, bty = 'n')
dev.off()

# export
write.table(stat_to_export, paste0(output_dir, "/stat_values_LTR_RT.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

save.image(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))

# load(paste0(output_dir, "/annotation_overlaps_per_base_venn_all_RM.RData"))