#!/usr/bin/env Rscript
library(rtracklayer)
library(reshape2)
library(stringr)
library(fmsb)
library(parallel)
library(igraph)
library(ggvenn)
library(ChIPpeakAnno)

source("compare_annot_utils.R")
# standard annot for evaluation:

## INPUT DATA
annot_std_gr <- import("reference_genomes/maize_B73/v3_MIPS_Repeats_loci_adjusted_names.gff3")
# are there any conflicts, oveverlapping annotations?
p <- findOverlaps(annot_std_gr, annot_std_gr)
ovlp <- p[!from(p) == to(p)]
# not empty - the are overlaps!!, must be cleaned up
# remove conflicts
annot_std_gr <- gff_cleanup(annot_std_gr) # this takes ~ 10-15 min
check_for_overlaps(annot_std_gr)


dante <- import("reference_genomes/maize_B73/B73_RefGen_v3.fa_dante.gff3", format = "gff3")
output_dir <- "maize_B73_plots"

annot_test_str <- c(
  DANTE_LTR="reference_genomes/maize_B73/B73_RefGen_v3.fa_dante_ltr.gff3",
  Inpactor2 = "reference_genomes/maize_B73/B73_RefGen_v3.fa_inpactor2_c3/Inpactor2_predictions.bed",
  EDTA = "reference_genomes/maize_B73/B73_RefGen_v3.fa.mod.EDTA.raw/B73_RefGen_v3.fa.mod.LTR.intact.gff3"
)
# genome size from ref seq
library(Biostrings)
s <- readDNAStringSet("reference_genomes/maize_B73/B73_RefGen_v3.fa")
# remove string after space in seq names
names(s) <- gsub(" .+$", "", names(s))
SL <- seqlengths(s)
# conversion table - rice to rexdb names

# info about TE protein domains:
te_info <- read.table("reference_genomes/lineage_domain_order.csv", sep="\t", header=TRUE)
te_domain_info <- strsplit(te_info$Domains.order, " ")
names(te_domain_info) <- gsub("_gypsy", "/gypsy", gsub("_copia", "/copia", gsub("/", "|", te_info$Lineage, fixed = TRUE)))


# genome size:
gs <- sum(SL)
# whole genome as GRanges
genome_gr <- GRanges(seqnames = names(SL), ranges = IRanges(start = 1, end = SL))

# add seq sizes to all annots
seqlengths(annot_std_gr) <- SL[seqlevels(annot_std_gr)]
seqlengths(genome_gr) <- SL[seqlevels(genome_gr)]


# invert intervals in annot_std_gr and adjust names
annot_str_gr_inv <- GenomicRanges::setdiff(genome_gr, annot_std_gr, ignore.strand = TRUE)
annot_str_gr_inv$Name <- "no_annotation"
annot_std <-  c(annot_std_gr, annot_str_gr_inv)
annot_std$ori_name <- annot_std$Name

# add seq sizes to all annots
annot_test_str_gr_raw <- sapply(annot_test_str, import)
for (i in names(annot_test_str_gr_raw)){
  seqlengths(annot_test_str_gr_raw[[i]]) <- SL[seqlevels(annot_test_str_gr_raw[[i]])]
}

# correct also attributes - Name should contain classification
annot_test_str_gr_raw <- sapply(annot_test_str_gr_raw, correct_attributes)


# calculate venn diagram per base:
annot_test_str_gr_raw_reduce <- sapply(annot_test_str_gr_raw, reduce)
pdf(paste0(output_dir, "/annotation_overlaps_per_base_venn.pdf"), width = 13, height = 13)
res <- makeVennDiagram(annot_test_str_gr_raw_reduce, by='base', NameOfPeaks = names(annot_test_str_gr_raw_reduce))
dev.off()
# export counts
VennCountsDF <- as.data.frame.matrix(res$vennCounts)
VennCountsDF$Perc <- VennCountsDF$Count / sum(VennCountsDF$Count) * 100
write.table(VennCountsDF, paste0(output_dir, "/annotation_overlaps_per_base_venn_all.csv"), sep = "\t", quote = FALSE, row.names = FALSE)



annot_test_str_gr_intact <- sapply(annot_test_str_gr_raw, get_only_intact)
annot_test_str_gr_intact <- sapply(annot_test_str_gr_intact, remove_overlaping_elements)

# complete means that all intervals are covered by annotations - no gaps
# gaps are now intervals labeled as "no_annotation"

annot_intact_complete <- sapply(annot_test_str_gr_intact, make_complete, genome_gr)
annot_intact_complete_adj <- sapply(annot_intact_complete, adjust_names)


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



# export annot pairs as csv table
write.table(annot_pairs$DANTE_LTR, paste0(output_dir, "/annot_pairs_DANTE_LTR.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(annot_pairs$Inpactor2, paste0(output_dir, "/annot_pairs_Inpactor2.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(annot_pairs$EDTA, paste0(output_dir, "/annot_pairs_EDTA.csv"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(stat_to_export, paste0(output_dir, "/stat_LTR_RT.csv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# data.frames for radarcharts:
# extract Class_I|LTR|Ty1/copia and Class_I|LTR|Ty3/gypsy
copia <- as.data.frame.matrix(t(sapply(annot_stat, function(x) unlist(x["Class_I|LTR|Ty1/copia",]))))[,5:10]
copia_df_export <- as.data.frame.matrix(t(sapply(annot_stat, function(x) unlist(x["Class_I|LTR|Ty1/copia",]))))
copia_df_export$Method <- rownames(copia_df_export)
copia_rc <-  rbind(rep(1, 6), rep(0,6), copia) # this is to set the scale
gypsy <- as.data.frame.matrix(t(sapply(annot_stat, function(x) unlist(x["Class_I|LTR|Ty3/gypsy",]))))[,5:10]
gypsy_df_export <- as.data.frame.matrix(t(sapply(annot_stat, function(x) unlist(x["Class_I|LTR|Ty3/gypsy",]))))
gypsy_df_export$Method <- rownames(gypsy_df_export)
gypsy_rc <-  rbind(rep(1, 6), rep(0,6), gypsy)

# compare elements using venn diagrams
annot_test_str_gr_copia <- sapply(annot_test_str_gr_intact, clean_gff, cls = "LTR/Copia")
annot_test_str_gr_gypsy <- sapply(annot_test_str_gr_intact, clean_gff, cls = "LTR/Gypsy")
annot_test_str_gr_ltr <- mapply(c, annot_test_str_gr_copia, annot_test_str_gr_gypsy, SIMPLIFY = FALSE)


copia_comparison <- compare_ltr(annot_test_str_gr_copia)
copia_annot_info <- append_dante_classification(copia_comparison$gr_all_unique, dante, te_domain_info)
copia_domains_completeness <- domain_completeness(copia_annot_info, cols = names(annot_test_str))
copia_comparison$venn_list_no_conflict <-  get_venn_without_domain_conflict(copia_annot_info,cols = names(annot_test_str))
copia_dante_classification <- class_summary(copia_annot_info, cols = names(annot_test_str))

gypsy_comparison <- compare_ltr(annot_test_str_gr_gypsy)
gypsy_annot_info <- append_dante_classification(gypsy_comparison$gr_all_unique, dante, te_domain_info)
gypsy_domains_completeness <- domain_completeness(gypsy_annot_info, cols = names(annot_test_str))
gypsy_comparison$venn_list_no_conflict <-  get_venn_without_domain_conflict(gypsy_annot_info,cols = names(annot_test_str))
gypsy_dante_classification <- class_summary(gypsy_annot_info, cols = names(annot_test_str))

copia_venn_count <- venn_counts_from_list(copia_comparison$venn_list)
gypsy_venn_count <- venn_counts_from_list(gypsy_comparison$venn_list)
gypsy_venn_count$domains <- factor("All", levels = c("No_conflict", "All"))
copia_venn_count$domains <- factor("All", levels = c("No_conflict", "All"))
copia_venn_count_good <- venn_counts_from_list(copia_comparison$venn_list_no_conflict)
gypsy_venn_count_good <- venn_counts_from_list(gypsy_comparison$venn_list_no_conflict)
copia_venn_count_good$domains <- factor("No_conflict", levels = c("No_conflict", "All"))
gypsy_venn_count_good$domains <- factor("No_conflict", levels = c("No_conflict", "All"))

copia_count_fin <- rbind(copia_venn_count,copia_venn_count_good)
gypsy_count_fin <- rbind(gypsy_venn_count, gypsy_venn_count_good)

# calcate the same but merge gypos and copia
ltr_comparison <- compare_ltr(annot_test_str_gr_ltr)
ltr_overlap_proportion <- get_overlap_proportion(annot_test_str_gr_ltr)
ltr_overlap_pass_threshold <- table(ltr_overlap_proportion >= 0.9)
ltr_annot_info <- append_dante_classification(ltr_comparison$gr_all_unique, dante, te_domain_info)
ltr_venn_count <- venn_counts_from_list(ltr_comparison$venn_list)
ltr_venn_count_plot <- ltr_venn_count[c("All", "DANTE_LTR", "EDTA", "Inpactor2"),]

# for main figure - LTR-RT elements (include oboth copia and gypsy)
pdf(paste0(output_dir, "/annotation_overlaps_ltr_rt.pdf"), width = 10, height = 5)
p1 <- ggvenn(ltr_comparison$venn_list,
             fill_color = c("#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF"), text_size=5) +
  theme(title = element_text(size = 30))
p2 <- ggplot(ltr_venn_count_plot, aes(x = method, y = counts)) +
  geom_bar(stat = "identity", position = "identity", fill = "#CCCCCC", color = "#000000") +
  labs(y = "Number of elements", x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    title = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.y = element_text(size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(linewidth = 2),
    plot.margin = unit(c(0.3, 1, 0.0, 1.5), "cm")
  )
multiplot(p2, p1, layout = matrix(c(0, 1,1,2, 2, 2), nrow=1))
dev.off()

gff2bases <- function (gr){
  grr <- reduce(gr)
  # create vector with seqid_base
  nt_index <- mapply(":", start(grr), end(grr))
  seq_nt <- mapply(paste0, as.vector(seqnames(grr)), nt_index)
}


# third panel
pdf(paste0(output_dir, "/annot_stat_comparison_full_elements_radarchart_ltr.pdf"), width = 5, height = 5, pointsize = 3)
colors <- c(DANTE_LTR = "red", Inpactor2 = "blue", EDTA = "green")
radarchart(stat_values_LTR_RT_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=3,
           title="", plwd = 4, pcol = c("red", "blue", "green", "black", "orange", "purple"))
# do not draw
legend(-1.4,-0.8, legend = names(colors), col = colors, lwd = 3, cex = 1.5, bty = 'n')
dev.off()


save.image("compare_annotation_maize1.RData")
# load("compare_annotation_maize1.RData")
# add also line labels

dir.create(output_dir, showWarnings = FALSE)
pdf(paste0(output_dir, "/annot_stat_comparison_full_elements_radarchart.pdf"),
    width = 16, height = 8, pointsize = 10)
par(mfrow = c(1,2), cex = 1, cex.main=2, mar = c(4,2,2,2))
colors <- c(DANTE_LTR = "red", Inpactor2 = "blue", EDTA = "green")
radarchart(copia_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=1.5,
           title="Ty1/copia", plwd = 3, pcol = colors[rownames(copia_rc)[-(1:2)]])
legend(-1,-1, legend = names(colors), col = colors, lwd = 3, cex = 1)
radarchart(gypsy_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=1.5,
           title="Ty3/gypsy", plwd = 3, pcol = colors[rownames(copia_rc)[-(1:2)]])
dev.off()


pdf(paste0(output_dir, "/annotation_overlaps.pdf"), width = 13, height = 13)
p1 <-  ggvenn(copia_comparison$venn_list) + ggtitle("Ty1/Copia") +
  theme(title = element_text(size = 20))
p2 <- ggvenn(copia_comparison$venn_list_no_conflict) + ggtitle("Ty1/Copia No conflict") +
  theme(title = element_text(size = 20))
p3 <- ggplot(copia_count_fin, aes(x = method, y = counts, fill=domains)) +
  geom_bar(stat = "identity", position = "identity", alpha=0.9, ) +
  # set the order of the fill levels
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(p1, p2, p3, cols = 2)

p1 <-  ggvenn(gypsy_comparison$venn_list) + ggtitle("Ty3/gypsy") +
  theme(title = element_text(size = 20))

p2 <- ggvenn(gypsy_comparison$venn_list_no_conflict) + ggtitle("Ty3/gypsy No conflict") +
    theme(title = element_text(size = 20))
p3 <- ggplot(gypsy_count_fin, aes(x = method, y = counts, fill=domains)) +
  geom_bar(stat = "identity", position = "identity", alpha=0.9, ) +
  # set the order of the fill levels
  theme(title = element_text(size = 10)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

multiplot(p1, p2, p3, cols = 2)
# set following plots p1 and p2 on the same page using ggarrange
dev.off()



# export tables - three groups:
write.table(copia_domains_completeness, paste0(output_dir, "/copia_domains_completeness.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(copia_dante_classification, paste0(output_dir, "/copia_dante_classification.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(gypsy_domains_completeness, paste0(output_dir, "/gypsy_domains_completeness.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(gypsy_dante_classification, paste0(output_dir, "/gypsy_dante_classification.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(copia_df_export, paste0(output_dir, "/copia_statistics.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(gypsy_df_export, paste0(output_dir, "/gypsy_statistics.csv"), sep="\t", quote=FALSE, row.names=FALSE)

# export annot_pairs as csv
write.table(annot_pairs$DANTE_LTR, paste0(output_dir, "/annot_pairs_DANTE_LTR.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(annot_pairs$Inpactor2, paste0(output_dir, "/annot_pairs_Inpactor2.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(annot_pairs$EDTA, paste0(output_dir, "/annot_pairs_EDTA.csv"), sep="\t", quote=FALSE, row.names=FALSE)