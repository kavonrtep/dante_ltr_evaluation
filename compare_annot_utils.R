multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

trim_gr <- function (grA, grMask){
  # this function trims genomicragnes in grA by genomicranges in grMask
  # metadata is preserved, if the resulting range is empty, reange is removed completelly
  # uses bedtools substract

  #grA = GRanges("chr1", IRanges(start = c(1, 10, 20, 30, 50, 60), width = c(5, 5, 5, 5, 5, 5)))
  #grA$Name= c("A", "B", "C", "D", "E", "F")
  #grMask = GRanges("chr1", IRanges(start = c(4, 20,55, 59), width = c(8,8,7,8)))

  tmp_gr <- tempfile(fileext = ".gff3")
  tmp_mask <- tempfile(fileext = ".gff3")
  tmp_out <- tempfile(fileext = ".gff3")

  export(grA, tmp_gr, format = "gff3")
  export(grMask, tmp_mask, format = "gff3")

  cmd <- paste0("bedtools subtract -a ", tmp_gr, " -b ", tmp_mask, " > ", tmp_out)
  system(cmd)

  gr_out <- import(tmp_out, format = "gff3")
  # delete tmp files
  unlink(tmp_gr)
  unlink(tmp_mask)
  unlink(tmp_out)
  return(gr_out)
}
clean_gff <- function(gr, cls){
  if (is.null(gr$source)){
    print(2)
    ## this is inpactor
    gr$name <- gsub("RLC", "LTR/Copia", gr$name)
    gr$name <- gsub("RLG", "LTR/Gypsy", gr$name)
    gr <- gr[startsWith(gr$name, cls)]
    #gr = reduce(gr)
    return(gr)
  }
  if (gr$source[1] == "EDTA"){
    print(1)
    gr <- gr[startsWith(gr$Classification, cls)]
    #gr = reduce(gr)
    return(gr)
  }

  if (any("dante_ltr" %in% gr$source)){
    print(3)
    gr$Final_Classification <- gsub("Class_I|LTR|Ty3/g", "LTR/Gypsy", gr$Final_Classification, fixed = TRUE)
    gr$Final_Classification <- gsub("Class_I|LTR|Ty1/c", "LTR/Copia", gr$Final_Classification, fixed = TRUE)
    gr$Final_Classification[is.na(gr$Final_Classification)] <- ""
    gr <- gr[startsWith(gr$Final_Classification, cls)]
    #gr = reduce(gr)
    return(gr)
  }
  if (gr$source[1] == 'RM'){
    print(5)
    gr$Name <- gsub("Class_I|LTR|Ty3/g", "LTR/Gypsy", gr$Name, fixed = TRUE)
    gr$Name <- gsub("Class_I|LTR|Ty1/c", "LTR/Copia", gr$Name, fixed = TRUE)
    gr$Name[is.na(gr$Name)] <- ""
    gr <- gr[startsWith(gr$Name, cls)]
    #gr = reduce(gr)
    return(gr)
  }
  print(4)
  gr <- gr[startsWith(gr$Name, cls)]
  #gr = reduce(gr)
  return(gr)
}
get_stat <- function(test, ref, genome_size){
  TP <- sum(width(GenomicRanges::intersect(test, ref, ignore.strand=TRUE)))
  FN <- sum(width(trim_gr(ref, test)))
  FP <- sum(width(trim_gr(test, ref)))
  TN <- genome_size - (TP + FN + FP)
  sens <- TP/(TP+FN)
  spec <- TN/(FP+TN)
  accu <- (TP+TN)/(TP+TN+FP+FN)
  prec <- TP/(TP+FP)
  FDR <- 1-prec
  F1 <- (2*TP)/(2 * TP+FP+FN)
  return(c(TP=TP,FN=FN, FP=FP, TN=TN, sens=sens, spec=spec,accu=accu,prec=prec, FDR=FDR, F1=F1))
}


get_FP <- function(test, ref){
  gr_fp <- trim_gr(test, ref)
  return(gr_fp)
}

plot_stats <- function(x, logax=''){
  y <- do.call(rbind, x)
  par(mfrow = c(2,1))
  barplot(y[,1:4], beside = TRUE, col=1:nrow(y), log=logax)
  legend('topleft', legend=rownames(y), col = 1:nrow(y),pch = 15, cex = 0.5)
  barplot(y[,-(1:4)], beside = TRUE, col = 1:nrow(y), ylim=c(0,1))
  abline(h=seq(0,1,0.1), lty=2)
}

get_only_intact <- function(gr){
  # keep only feature for full elements, remove subfeatures
  if (!is.null(gr$Rank)){
    return(gr[gr$Rank != "D" & gr$type=="transposable_element"])
  }
  if (is.null(gr$source)){
    # inpactor
    return(gr)
  }
  if (gr$source[1]=="EDTA"){
   return(gr[gr$type=="repeat_region"])
  }
}


## compare_ltr = function(gr_list){
##   gr_all = reduce(unlist(GRangesList(gr_list)))
##   ovlp = matrix(data=0, ncol=length(gr_list), nrow=length(gr_all))
##   colnames(ovlp) = names(gr_list)
##   for (i in names(gr_list)){
##     ovlp[to(findOverlaps(gr_list[[i]], gr_all, minoverlap = 1500)),i] = 1
##   }
##   ovlp_domains = rep(0, length(gr_all))
##                                         # count domains
##   Nd = tabulate(to(findOverlaps(dante, gr_all)))
##   ovlp_domains[1:length(Nd)] = Nd
## }

compare_two_ltr_sets <- function(gr1, gr2){
  p <- findOverlaps(gr1, gr2)
  w1 <- width(gr1[from(p)])
  w2 <- width(gr2[to(p)])
  p_ovl_w <- findOverlapPairs(gr1, gr2)
  ovl_perc <- width(pintersect(p_ovl_w))/ifelse(w1>w2, w1, w2)
  N12 <- sum(ovl_perc>0.95)
  N1 <- length(gr1) - N12
  N2 <- length(gr2) - N12
  return(c(N1=N1,N2=N2,N12=N12))
}

get_ranges_from_cluster <- function(x){
  p <- findOverlapPairs(x)
  pi <- findOverlaps(x)
  W1 <- width(x[from(pi)])
  W2 <- width(x[to(pi)])
  Wmin <- ifelse(W1>W2, W2, W1)
  W <- width(pintersect(p))
  to_merge <- pi[W/Wmin > 0.9]
  cls <- clusters(graph_from_data_frame(as.data.frame(to_merge)))$membership
  gout <- unlist(GRangesList(sapply(split(x, cls), GenomicRanges::reduce, ignore.strand=TRUE)))
}

get_overlap_proportion <- function(gr_list){
  grc <- unlist(GRangesList(gr_list))
  ovlp <- findOverlapPairs(grc)
  s1 <- start(ovlp@first)
  e1 <- end(ovlp@first)
  s2 <- start(ovlp@second)
  e2 <- end(ovlp@second)
  # calculate overlap proportion os percent of overlap size over shorted interval, use start and end positions
  ovlp_prop <- width(pintersect(ovlp))/ifelse((e1-s1) < (e2 - s2), e1-s1 + 1, e2-s2 + 1)
  ovlp_prop
}

compare_ltr <- function(gr_list){
  grc <- unlist(GRangesList(gr_list))
  gr_all <- GenomicRanges::reduce(grc, with.revmap=TRUE, ignore.strand=TRUE)
  ## analyze clusters:
  singletons <- gr_all[sapply(gr_all$revmap, length)==1]
  gr_cls <- gr_all[sapply(gr_all$revmap, length)>1]
  cls <- gr_all$revmap[sapply(gr_all$revmap, length)>1]
  ## fully or nearly fully overlapping
  ## good_to_merge = sapply(cls, function(x) min(width(grc[x])))/width(gr_cls)>0.9
  W <- width(gr_cls[seq_along(cls)])
  # merge is overlap is at least 90% of cluster size
  good_to_merge <- sapply(seq_along(cls), function(x) min(width(grc[cls[[x]]]))/W[x]>0.9)
  gr_out1 <- gr_cls[good_to_merge]
  cluster_to_inspect <- gr_cls[!good_to_merge]

  ## all overlaps:
  gr_out2 <- unlist(GRangesList(sapply(sapply(cluster_to_inspect$revmap, function(x)grc[x]), get_ranges_from_cluster)))
  gr_all_unique <- c(singletons, gr_out1, gr_out2)
  ## get overlaps
  hit_all <- list()
  for (i in seq_along(gr_list)){
    p <- findOverlaps(gr_all_unique, gr_list[[i]], ignore.strand=TRUE)
    w1 <- width(gr_all_unique[from(p)])
    w2 <- width(gr_list[[i]][to(p)])
    p_ovl_w <- findOverlapPairs(gr_all_unique, gr_list[[i]])
    ovl_perc <- width(pintersect(p_ovl_w))/ifelse(w1>w2, w1, w2)
    hit <- numeric(length(gr_all_unique))
    hit[from(p)[ovl_perc>0.9]] <- 1
    hit_all[[names(gr_list[i])]] <- hit
    mcols(gr_all_unique)[names(gr_list)[i]] <- hit
  }
  hits <- do.call(cbind, hit_all)
  venn_list <- apply(hits, 2, function(x)which(x==1))
  return(list(gr_all_unique = gr_all_unique, hits = hits, venn_list=venn_list))
}

resolve_name <- function(x){
  if (length(x)==1){
                                        # no conflict
    return(x)
  } else{
    y <- sapply(x, strsplit, split="|", fixed = TRUE)
    ny <- table(unlist(sapply(y, function(x)paste(seq_along(x), x))))
    if (max(ny)<length(x)){
      return("Unknown")
    }else{
      k <- which(ny==length(x))
      r <- max(as.numeric((gsub(" .+", "", names(k)))))
      out <- paste(y[[1]][1:r], collapse="|")
      return(out)
    }
  }
}

resolve_name_wicker <- function(x){
  if (length(x)==1){
    return(x)
  } else{
    y <- sapply(x, strsplit, split="", fixed = TRUE)
    ny <- table(unlist(sapply(y, function(x)paste(seq_along(x), x))))
    if (max(ny)<length(x)){
      return("XXX")
    }else{
      k <- which(ny==length(x))
      r <- max(as.numeric((gsub(" .+", "", names(k)))))
      out <- paste(y[[1]][1:r], collapse="")
      if (out == "-"){
        return("-")
      }
      if (nchar(out)<3) {
        out <- paste0(out, paste(rep("X", 3-nchar(out)), collapse=""))
      }
      return(out)
    }
  }
}

dante_filtering <- function(dante_gff, min_similarity=0.4,
                            min_identity=0.2, Relative_Length=0.6,
                            min_relat_interuptions=8) {
  include <- as.numeric(dante_gff$Similarity) >= min_similarity &
    as.numeric(dante_gff$Identity) >= min_identity &
    as.numeric(dante_gff$Relat_Length) >= Relative_Length &
    as.numeric(dante_gff$Relat_Interruptions) <= min_relat_interuptions
  include[is.na(include)] <- FALSE

  return(dante_gff[include])
}


append_dante_classification <- function(gr, dante, te_domain_info){
  ## filter dante --
  dante_f <- dante_filtering(dante)
  #ovl = findOverlaps(gr, dante_f, type = 'within')
  ovl <- findOverlaps(dante_f, gr, type ='within' )
  #dante_groups = split(dante_f$Final_Classification[to(ovl)], from(ovl))
  dante_groups <- split(dante_f$Final_Classification[from(ovl)], to(ovl))
  domains_groups <- split(dante_f$Name[from(ovl)], to(ovl))
  lca_annot <- sapply(lapply(dante_groups, unique), resolve_name)
  domain_completenes <- mapply(function(x, y){sum(y %in% te_domain_info[[x]])/length(te_domain_info[[x]])}, lca_annot, domains_groups)
  domain_completenes[is.nan(domain_completenes)] <- 0
  domain_duplication <- sapply(domains_groups, function(x)sum(duplicated(x)))
  domain_names <- sapply(domains_groups, function(x)paste(x, collapse = ","))
  ndomains <- sapply(dante_groups, length)
  gr$classification <- "no_domain"
  gr$ndomains <- 0
  gr$domain_completeness <- 0
  gr$domain_duplication <- 0
  gr$domain_names <- ""
  gr$classification[as.numeric(names(lca_annot))] <- lca_annot
  gr$ndomains[as.numeric(names(ndomains))] <- ndomains
  gr$domain_completeness[as.numeric(names(ndomains))] <- domain_completenes
  gr$domain_duplication[as.numeric(names(ndomains))] <- domain_duplication
  gr$domain_names[as.numeric(names(ndomains))] <- domain_names

  return(gr)

}

venn_counts_from_list <- function(venn_list){
  N <- length(venn_list)
  counts <- c(All = length(unique(unlist(venn_list))))
  for (i in 1:(N-1)){
    label <- combn(names(venn_list), i)
    for (j in 1:ncol(label)){
      counts[paste(label[,j], collapse = ",")] <- length(unique(unlist(venn_list[label[,j]])))
    }
  }
  counts <-  data.frame(counts=counts)
  counts$method <-  rownames(counts)
  return(counts)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

class_summary <- function(info, cols = c('edta', 'dante_ltr', 'dante_ltr_clean', 'inpactor')){
  cls <- sort(unique(info$classification))
  cls_counts <- c(table(factor(info$classification, levels = cls)))
  df <- data.frame(Annotation=cls, total = cls_counts)
  for (n in cols){
    df[[n]] <- c(table(factor(info$classification[mcols(info)[, n]==1], levels = cls)))
  }
  return(df)
}

domain_completeness <- function(info, cols){
  dc <- list()
  for (n in cols){
    x <- info$domain_completeness[mcols(info)[, n] == 1]
    x[is.na(x)] <- 0
    conflict <- sum(x == 0 & info$ndomains[mcols(info)[, n] == 1]>0)  # domain completeness is zero because of conflict classification
    dc[[n]] <- c(
      conflict = conflict,
      absent = sum(x==0) - conflict,
      partially_missing = sum(x>0 & x<1),
      complete = sum(x>=1 & x < 1.17),
      duplicated = sum(x >= 1.17)
    )
  }
  out <- as.data.frame(do.call(rbind, dc))
  out <- cbind(method = names(dc), out)
  return(out)
}

get_venn_without_domain_conflict <- function(info, cols){
  is_conflict <- info$domain_completeness == 0 & info$ndomains > 0
  vl <- list()
  for (n in cols){
    vl[[n]] <- which(mcols(info)[, n] == 1 & !is_conflict)
  }
  return(vl)
}

adjust_names <- function (g){
  if (is.null(g$source)){
    #  it is Inpactor2
    new_name <- ifelse(grepl("RLC", g$Name), "Class_I|LTR|Ty1/copia", g$Name)
    new_name <- ifelse(grepl("RLG", new_name), "Class_I|LTR|Ty3/gypsy", new_name)
    g$Name <- new_name
    print("inpactor")
    return(g)
  }
  if (g$source[1] == "Inpactor2"){
    new_name <- ifelse(grepl("RLC", g$Name), "Class_I|LTR|Ty1/copia", g$Name)
    new_name <- ifelse(grepl("RLG", new_name), "Class_I|LTR|Ty3/gypsy", new_name)
    g$Name <- new_name
    print("inpactor2")
    return(g)
  }
  if (g$source[1] == 'dante_ltr'){
    # truncate dante names
    g$ori_name <- g$Name
    new_name <- gsub("copia[|].*", "copia", g$Name)
    new_name <- gsub("gypsy[|].*", "gypsy", new_name)
    g$Name <- new_name
    print("dante_ltr")
    return(g)
  }
  if (g$source[1] == 'EDTA'){
    # truncate dante names
    new_name <- ifelse(g$Name == 'LTR/Gypsy', 'Class_I|LTR|Ty3/gypsy', g$Name)
    new_name <- ifelse(new_name == 'LTR/Copia', 'Class_I|LTR|Ty1/copia', new_name)
    new_name <- ifelse(new_name == 'LTR/unknown', 'Class_I|LTR', new_name)
    g$Name <- new_name
    print("EDTA")
    return(g)
  }
  print("unknown source")
}

adjust_names2 <- function (g){
  if (is.null(g$source)){
    #  it is Inpactor2
    new_name <- ifelse(grepl("RLC", g$Name), "Class_I|LTR|Ty1/copia", g$Name)
    new_name <- ifelse(grepl("RLG", new_name), "Class_I|LTR|Ty3/gypsy", new_name)
    g$Name <- new_name
    print("inpactor")
    return(g)
  }
  if (g$source[1] == "Inpactor2"){
    new_name <- ifelse(grepl("RLC", g$Name), "Class_I|LTR|Ty1/copia", g$Name)
    new_name <- ifelse(grepl("RLG", new_name), "Class_I|LTR|Ty3/gypsy", new_name)
    g$Name <- new_name
    print("inpactor2")
    return(g)
  }
  if (g$source[1] == 'dante_ltr'){
    # truncate dante names
    g$ori_name <- g$Name
    N <- gsub("-","/", gsub("[/]", "|", g$Name))
    new_name <- gsub("copia[|].*", "copia", N)
    new_name <- gsub("gypsy[|].*", "gypsy", new_name)
    g$Name <- new_name
    print("dante_ltr")
    return(g)
  }
  if (g$source[1] == 'EDTA'){
    # truncate dante names
    new_name <- ifelse(g$Name == 'LTR/Gypsy', 'Class_I|LTR|Ty3/gypsy', g$Name)
    new_name <- ifelse(new_name == 'LTR/Copia', 'Class_I|LTR|Ty1/copia', new_name)
    new_name <- ifelse(new_name == 'LTR/unknown', 'Class_I|LTR', new_name)
    g$Name <- new_name
    print("EDTA")
    return(g)
  }
  print("unknown source")
}


correct_attributes <- function(g){
  if (is.null(g$source)){
    #  it is Inpactor2
    g$Name <- g$name
    return(g)
  }
  if (g$source[1] == 'dante_ltr'){
    # truncate dante names
    return(g)
  }
  if (g$source[1] == 'EDTA'){
    # truncate dante names
    g$Name <- g$Classification
    return(g)
  }
}



gff_cleanup <- function(gff, wicker = FALSE){
  ## remove overlapin annotation track - assign new annot
  if (wicker){
    resolve_name_fun <- resolve_name_wicker
  }else{
    resolve_name_fun <- resolve_name
  }
  gff_disjoin <- disjoin(gff, with.revmap=TRUE, ignore.strand=TRUE)
  ## append annotation:
  gff_names <- mclapply(as.list(gff_disjoin$revmap), FUN = function(x)gff$Name[x], mc.cores = 8)
  gff_strands <- mclapply(as.list(gff_disjoin$revmap), FUN = function(x)strand(gff[x]), mc.cores = 8)
  new_annot <- sapply(sapply(gff_names, unique), paste, collapse=":")
  new_annot_uniq <- unique(new_annot)
  lca_annot <- sapply(strsplit(new_annot_uniq, ":", fixed = TRUE), resolve_name_fun)
  names(lca_annot) <- new_annot_uniq
  new_annot_lca <- lca_annot[new_annot]
  #new_annot_lca = sapply(sapply(gff_names, unique), resolve_name)
  strand_attribute <- sapply(sapply(gff_strands, unique), paste, collapse="|")
  gff_disjoin$source <- "RM"
  gff_disjoin$type <- "repeat"
  gff_disjoin$score <- NA
  gff_disjoin$phase <- NA
  gff_disjoin$Name <- new_annot_lca
  gff_disjoin$Original_names <- new_annot
  gff_disjoin$strands <- strand_attribute
  gff_disjoin$revmap <- NULL
  return(gff_disjoin)
}

check_for_overlaps <- function(g){
  # check for overlaps
  ovl <- findOverlaps(g, g, ignore.strand = TRUE)
  ovl <- ovl[from(ovl) != to(ovl)]
  if (length(ovl) > 0){
    print("overlaps found")
    return(TRUE)
  } else {
    print("no overlaps found")
    return(FALSE)
  }
}

remove_overlaping_elements <- function(g){
  # remove overlaping elements
  ovl <- findOverlaps(g, g, ignore.strand = TRUE)
  ovl <- ovl[from(ovl) != to(ovl)]
  to_rm <- unique(to(ovl))
  if (length(to_rm) > 0){
    print(paste0("removing ", length(to_rm), " elements"))
    g <- g[-to_rm]
  }
    g
}

get_overlaping_elements <- function(g){
    # remove overlaping elements
    ovl <- findOverlaps(g, g, ignore.strand = TRUE)
    ovl <- ovl[from(ovl) != to(ovl)]
    to_rm <- unique(to(ovl))
    if (length(to_rm) > 0){
        print(paste0("removing ", length(to_rm), " elements"))
        g <- g[to_rm]
    }
    g

}


make_complete <- function(g, genome_gr){
  # fill in missing intervals, i.e. intervals without annotation
  g_m <- GenomicRanges::setdiff(genome_gr, g, ignore.strand = TRUE)
  g_m$Name <- "no_annotation"
  g_c <- append(g, g_m)
  g_c
}

make_complete2 <- function(g, s){
  SL <- seqlengths(s)
  genome_gr <- GRanges(seqnames = names(SL), ranges = IRanges(start = 1, end = SL))
    # fill in missing intervals, i.e. intervals without annotation
  g_m <- GenomicRanges::setdiff(genome_gr, g, ignore.strand = TRUE)
  g_m$Name <- "no_annotation"
  g_c <- append(g, g_m)
}




extract_all_matches_as_gr <- function(annot_test, annot_ref){
  g12 <- append(annot_test, annot_ref)
  g12_dis <-  GenomicRanges::disjoin(g12, ignore.strand = TRUE)
  # for each element in g12_dis find the corresponding element in g1
  g1ovl <- findOverlaps(g12_dis, annot_ref, ignore.strand = TRUE)
  g2ovl <- findOverlaps(g12_dis, annot_test, ignore.strand = TRUE)
  # check if all is matched exactly once
  if (all(from(g1ovl) == from(g2ovl))){
    print("all elements matched exactly once")
  }else{
    print("not all elements matched exactly once")
  }
  g12_dis$ref <- annot_ref[to(g1ovl)]$Name
  g12_dis$test <- annot_test[to(g2ovl)]$Name
  g12_dis$Name <- paste0(g12_dis$ref, ":", g12_dis$test)
  return(g12_dis)
}
sanitize_filename <- function(filename) {
  # Replace special characters with a safe alternative
  sanitized <- gsub("[/\\\\:*?\"<>|]", "-", filename)
  return(sanitized)
}


extract_all_matches <- function(annot_test, annot_ref){
  # annot_ref is the reference
  # annot_test is the test
  # output is a data.frame with all identified matching pairs and the number of basepairs
  g12 <- append(annot_test, annot_ref)
  g12_dis <-  GenomicRanges::disjoin(g12, ignore.strand = TRUE)
  # for each element in g12_dis find the corresponding element in g1
  g1ovl <- findOverlaps(g12_dis, annot_ref, ignore.strand = TRUE)
  g2ovl <- findOverlaps(g12_dis, annot_test, ignore.strand = TRUE)
  # check if all is matched exactly once
  if (all(from(g1ovl) == from(g2ovl))){
    print("all elements matched exactly once")
  }else{
    print("not all elements matched exactly once")
  }
  N1 <- annot_ref[to(g1ovl)]$Name
  N2 <- annot_test[to(g2ovl)]$Name
  N1_lengths <- width(g12_dis[from(g1ovl)])
  N2_lengths <- width(g12_dis[from(g2ovl)])
  tabulated_data <- aggregate(N1_lengths, by = list(N1, N2), FUN = sum)
  colnames(tabulated_data) <- c('ref', 'test', 'basepairs')
  return(tabulated_data)
}

get_matched_granges <- function(annot_test, annot_ref){
  g12 <- append(annot_test, annot_ref)
  g12_dis <-  GenomicRanges::disjoin(g12, ignore.strand = TRUE)
  # for each element in g12_dis find the corresponding element in g1
  g1ovl <- findOverlaps(g12_dis, annot_ref, ignore.strand = TRUE)
  g2ovl <- findOverlaps(g12_dis, annot_test, ignore.strand = TRUE)
  g12_dis$ref <- annot_ref[to(g1ovl)]$Name
  g12_dis$test <- annot_test[to(g2ovl)]$Name
  g12_dis$Name <- paste0(g12_dis$test, ":", g12_dis$ref)
  split_gr <- split(g12_dis, g12_dis$Name)
  return(split_gr)
}


evaluate_matches <- function (annot_pairs){
  Names <- annot_pairs[, 1]
  matches <- str_detect(annot_pairs[,2], fixed(Names))
  # NA in matches means that there was conflict in ref annot.
  # true positives. must match, and not be no_annotation and not be NA
  TP <- matches; TP[is.na(TP)] <- FALSE; TP[Names == "no_annotation"] <- FALSE
  # true negatives
  TN <- annot_pairs[,1] == "no_annotation" & annot_pairs[,2] == "no_annotation"
  # false negative -
  FN <- Names != 'no_annotation' & !matches; FN[is.na(FN)] <- FALSE
  # false positive
  FP <- Names == 'no_annotation' & !matches; FP[is.na(FP)] <- FALSE
  anno_pairs <- data.frame(annot_pairs, matches, TP, TN, FN, FP)
  lvl <- sort(unique(c(annot_pairs[, 1], annot_pairs[, 2])))
  anno_pairs$ref <- factor(anno_pairs$ref, levels = lvl)
  anno_pairs$test <- factor(anno_pairs$test, levels = lvl)
  return(anno_pairs)
}

get_confusion_matrix <- function(annot_pairs){
  p1 <- evaluate_matches(annot_pairs)
  p1m <- dcast(p1, ref ~ test, value.var = "basepairs", drop = FALSE)
  rownames(p1m) <- p1m$ref
  conf_matrix <-  (as.matrix(p1m[,-1]))
  conf_matrix[is.na(conf_matrix)] <- 0
  return(conf_matrix)
}

get_stat_from_confusion_matrix <- function (cm){
    # cm is a confusion matrix
    # returns a data.frame with the following columns:
    # TP, FP, FN, TN, precision, recall, F1, accuracy
    TP <- diag(cm)
    FP <- colSums(cm) - TP
    FN <- rowSums(cm) - TP
    TN <- sum(cm) - TP - FP - FN
    precision <- TP / (TP + FP)
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    recall <- sensitivity

    F1 <- 2 * precision * recall / (precision + recall)
    accuracy <- (TP + TN) / (TP + TN + FP + FN)
    FDR <- FP / (TP + FP)
    stat <- data.frame(TP, FP, FN, TN, precision, sensitivity,specificity, FDR, F1, accuracy)
    return(stat)
}


calculate_TP_TF_FP_FN_from_pairs_LTR_category <- function(ap){
  tp_cases <- read.table("/mnt/raid/454_data/dante/reference_genomes/TP_FP_cases.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  error_type <- tp_cases$stat
  names(error_type) <- paste(tp_cases$ref, tp_cases$test, sep = "__")
  stat_call <- c(TP=0, FN=0, FP=0, TN=0)
  for (i in 1:nrow(ap)){
    pair <- paste(ap$ref[i], ap$test[i], sep = "__")
    e_type <- error_type[pair]
    if (is.na(e_type)){
      print(pair)
      stop("error type not found")

    }
    stat_call[e_type] <- stat_call[e_type] + ap$basepairs[i]
    print("xxxxxxxxxxxxxxxxx")
    print(i)
    print(e_type)
    print(stat_call)
    print("-----------------")
  }
  return(stat_call)
}

calculate_statistics_from_groups <- function(stat_call){
  # calculate precision sensitivity specificity        FDR        F1  accuracy
  # from TP, TF, FP, FN
  out <- c(
    precision = unname(stat_call["TP"] / (stat_call["TP"] + stat_call["FP"])),
    sensitivity = unname(stat_call["TP"] / (stat_call["TP"] + stat_call["FN"])),
    specificity = unname(stat_call["TN"] / (stat_call["TN"] + stat_call["FP"])),
    FDR = unname(stat_call["FP"] / (stat_call["FP"] + stat_call["TP"])),
    F1 = unname(2 * stat_call["TP"] / (2 * stat_call["TP"] + stat_call["FP"] + stat_call["FN"])),
    accuracy = unname((stat_call["TP"] + stat_call["TN"]) / sum(stat_call))
  )
  return(out)

}


