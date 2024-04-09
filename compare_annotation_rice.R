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
annot_std_gr <- import("rice_v7/riceTElib/rice_v7.fasta.gff3")
dante <- import("/mnt/raid/454_data/dante/reference_genomes/rice_v7/rice_v7.fasta_dante.gff3", format = "gff3")
output_dir <- "rice_v7_plots"

annot_test_str <- c(
  DANTE_LTR="/mnt/raid/454_data/dante/reference_genomes/rice_v7/rice_v7.fasta_dante_ltr.gff3",
  Inpactor2 = "/mnt/raid/454_data/dante/reference_genomes/rice_v7/rice_v7.fasta_inpactor2_c3/Inpactor2_predictions.bed",
  EDTA = "/mnt/raid/454_data/dante/reference_genomes/rice_v7/rice_v7.fasta.mod.EDTA.raw/rice_v7.fasta.mod.LTR.intact.gff3"
)
# genome size from ref seq
library(Biostrings)
s <- readDNAStringSet("rice_v7/rice_v7.fasta")
SL <- seqlengths(s)
# conversion table - rice to rexdb names
rice2rexdb <- read.table(
  "rice_v7/riceTElib//rice7.0.0.liban_unique_categories_to_rexdb.csv",
  sep = "\t", header = FALSE, as.is = TRUE)
rownames(rice2rexdb) <- rice2rexdb[,1]

# info about TE protein domains:
te_info <- read.table("/mnt/raid/users/petr/workspace/dante_ltr/databases/lineage_domain_order.csv", sep="\t", header=TRUE)
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
annot_str_gr_inv <- setdiff(genome_gr, annot_std_gr, ignore.strand = TRUE)
annot_str_gr_inv$Name <- "no_annotation"
annot_std <-  c(annot_std_gr, annot_str_gr_inv)
annot_std$ori_name <- annot_std$Name
annot_std$Name <- rice2rexdb[annot_std$Name,2]
annot_std$Name[is.na(annot_std$Name)] <- "Uknown"

# add seq sizes to all annots
annot_test_str_gr_raw <- sapply(annot_test_str, import)
for (i in names(annot_test_str_gr_raw)){
  seqlengths(annot_test_str_gr_raw[[i]]) <- SL[seqlevels(annot_test_str_gr_raw[[i]])]
}

# calculate venn diagram per base:
annot_test_str_gr_raw_reduce <- sapply(annot_test_str_gr_raw, reduce)
pdf(paste0(output_dir, "/annotation_overlaps_per_base_venn.pdf"), width = 13, height = 13)
res <- makeVennDiagram(annot_test_str_gr_raw_reduce, by='base', NameOfPeaks = names(annot_test_str_gr_raw_reduce))
dev.off()
# export counts
VennCountsDF <- as.data.frame.matrix(res$vennCounts)
VennCountsDF$Perc <- VennCountsDF$Count / sum(VennCountsDF$Count) * 100

write.table(VennCountsDF, paste0(output_dir, "/annotation_overlaps_per_base_venn_all.csv"), sep = "\t", quote = FALSE, row.names = FALSE)





# correct also attributes - Name should contain classification
annot_test_str_gr_raw <- sapply(annot_test_str_gr_raw, correct_attributes)


annot_test_str_gr_intact <- sapply(annot_test_str_gr_raw, get_only_intact)
annot_test_str_gr_intact <- sapply(annot_test_str_gr_intact, remove_overlaping_elements)
annot_intact_complete <- sapply(annot_test_str_gr_intact, make_complete, genome_gr) # this fill gaps between elements as no_annotation
annot_intact_complete_adj <- sapply(annot_intact_complete, adjust_names)
# adjust names
# only names matching the reference


annot_pairs <- lapply(annot_intact_complete_adj, extract_all_matches,annot_std)  # compare to reference annnotation
annot_confusion_matrix <- lapply(annot_pairs, get_confusion_matrix)
annot_stat <- lapply(annot_confusion_matrix, get_stat_from_confusion_matrix)

# comparison of all LTR-RT elements, for the purpose of comparison, copia annotated as gypsy is considered false positive
stat_LTR_RT <- lapply(annot_pairs, calculate_TP_TF_FP_FN_from_pairs_LTR_category)
stat_values_LTR_RT <- lapply(stat_LTR_RT, calculate_statistics_from_groups)
stat_values_LTR_RT <- as.data.frame.matrix(do.call(rbind, stat_values_LTR_RT))
stat_values_LTR_RT_rc <-  rbind(rep(1, 6), rep(0,6), stat_values_LTR_RT) # this is to set the scale
stat_to_export <- cbind(do.call(rbind, stat_LTR_RT), stat_values_LTR_RT)
stat_to_export$Method <- names(stat_LTR_RT)

# export annot pairs a csv table
write.table(annot_pairs$DANTE_LTR, file = "rice_v7_plots/annot_pairs_DANTE_LTR.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(annot_pairs$Inpactor2, file = "rice_v7_plots/annot_pairs_Inpactor2.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(annot_pairs$EDTA, file = "rice_v7_plots/annot_pairs_EDTA.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# export stat values for all LTR-RT elements
write.table(stat_to_export, file = "rice_v7_plots/stat_values_LTR_RT.csv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



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
annot_test_str_gr_copia <- sapply(annot_test_str_gr_intact, clean_gff, cls = "LTR/Copia") # this adjust attributes and names to be same in all annotations
annot_test_str_gr_gypsy <- sapply(annot_test_str_gr_intact, clean_gff, cls = "LTR/Gypsy")
annot_test_str_gr_ltr <- mapply(c, annot_test_str_gr_copia, annot_test_str_gr_gypsy, SIMPLIFY = FALSE)


# *compare_ltr* compare all elements in all annotations for venn diagrams
copia_comparison <- compare_ltr(annot_test_str_gr_copia)
copia_overlap_proportion <- get_overlap_proportion(annot_test_str_gr_copia)
copia_annot_info <- append_dante_classification(copia_comparison$gr_all_unique, dante, te_domain_info)
copia_domains_completeness <- domain_completeness(copia_annot_info, cols = names(annot_test_str))
copia_comparison$venn_list_no_conflict <-  get_venn_without_domain_conflict(copia_annot_info,cols = names(annot_test_str))
copia_dante_classification <- class_summary(copia_annot_info, cols = names(annot_test_str))

gypsy_comparison <- compare_ltr(annot_test_str_gr_gypsy)
gypsy_overlap_proportion <- get_overlap_proportion(annot_test_str_gr_gypsy)
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

# calculate the same but merge gypsy and copia
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





save.image("compare_annotation_rice1.RData")
# load("compare_annotation_rice1.RData")

# third panel
pdf("rice_v7_plots/annot_stat_comparison_full_elements_radarchart_ltr.pdf", width = 5, height = 5, pointsize = 3)
colors <- c(DANTE_LTR = "red", Inpactor2 = "blue", EDTA = "green")
radarchart(stat_values_LTR_RT_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=3,
           title="", plwd = 4, pcol = c("red", "blue", "green", "black", "orange", "purple"))
# do not draw
legend(-1.4,-0.8, legend = names(colors), col = colors, lwd = 3, cex = 1.5, bty = 'n')
dev.off()



## PLOTS
# add also line labels
# this is supplementary figure
dir.create(output_dir, showWarnings = FALSE)
pdf("rice_v7_plots/annot_stat_comparison_full_elements_radarchart_copia_gypsy.pdf",
    width = 16, height = 8, pointsize = 4)
par(mfrow = c(1,2), cex = 1, cex.main=5, mar = c(4,2,2,2))
colors <- c(DANTE_LTR = "red", Inpactor2 = "blue", EDTA = "green")
radarchart(copia_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=3,
           title="Ty1/copia", plwd = 4, pcol = colors[rownames(copia_rc)[-(1:2)]])
legend(-1,-1, legend = names(colors), col = colors, lwd = 3, cex = 2)
radarchart(gypsy_rc, axistype=1, seg=5, caxislabels = seq(0,1,.2), plty=1, vlcex=3,
           title="Ty3/gypsy", plwd = 4, pcol = colors[rownames(copia_rc)[-(1:2)]])
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


# export tables:
write.table(copia_domains_completeness, paste0(output_dir, "/copia_domains_completeness.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(copia_dante_classification, paste0(output_dir, "/copia_dante_classification.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(gypsy_domains_completeness, paste0(output_dir, "/gypsy_domains_completeness.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(gypsy_dante_classification, paste0(output_dir, "/gypsy_dante_classification.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(copia_df_export, paste0(output_dir, "/copia_statistics.csv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(gypsy_df_export, paste0(output_dir, "/gypsy_statistics.csv"), sep="\t", quote=FALSE, row.names=FALSE)
