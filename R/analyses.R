#' @title Intersection analysis
#' @description Calculate if overlap between different gene list ids is significant
#' @author Alexander Gabel
#' @usage intersect.analysis(gene.list, geneUniverse, alternative="greater")
#' @param gene.list list containing vectors of gene ids
#' @param geneUniverse all gene ids in the dataset
#' @param alternative character string indicating the alternative hypothesis and must be one of "two.sided", "greater" (default), or "less". Details \link{fisher.test}
#' @return a \code{plot} 
#' @export
intersect.analysis <- function(gene.list, geneUniverse, alternative="greater"){
  
  list.names <- names(gene.list)
  
  p.values <- c()
  p.names <- c()
  
  A <- B <- c()
  
  if(length(B) > length(A)){
    A <- sort(gene.list[[2]])
    B <- sort(gene.list[[1]])
  }else{
    A <- sort(gene.list[[1]])
    B <- sort(gene.list[[2]])
  }
  n <- length(geneUniverse) - length(which(!A %in% geneUniverse))- length(which(!B %in% geneUniverse))
  
  contingency.matrix <- matrix(c(n - length(union(A,B)), length(setdiff(A,B)), length(setdiff(B,A)), length(intersect(A,B))), nrow=2)
  fisher.results <- fisher.test(contingency.matrix, alternative=alternative)
  p.values <- fisher.results$p.value
  p.names <- paste0(list.names[1],"_vs_",list.names[2])
  
  names(p.values) <- p.names
  
  result <- c(length(intersect(A,B)), p.values)
  names(result) <- c(paste0(c("Count:", "P.value:"),c(p.names,p.names)))
  
  return(list(stat=result, genes=intersect(A,B)))
}

#' @title Fisher test for transcriptional overlap
#' @description Calculate significance if ratio of differentially expressed genes from a given gene set is significantly different to the ratio of differentially expressed genes in the complete dataset
#' @author Alexander Gabel
#' @usage do.fisher.test(fc.list, fg.ids, ml.list.up, ml.list.down, fc.threshold, ml.threshold=0.9, alternative = "two.sided")
#' @param fc.list list containing the log2-foldchange matrices of treatment [grafted|separated] [top|bottom] vs. intact  
#' @param fg.ids vector containing gene ids which shall be used for overlap analysis 
#' @param ml.list.up list containing matrices with marginal likelihoods if a gene is up regulated
#' @param ml.list.down list containing matrices with marginal likelihoods if a gene is down regulated
#' @param fc.threshold log2 foldchange threshold to define if a gene is differentially expressed
#' @param ml.threshold marginal likelihood threshold to define if a gene is differentially expressed
#' @param alternative character string indicating the alternative hypothesis and must be one of "two.sided" (default), "greater", or "less". Details \link{fisher.test}
#' @return a \code{list} 
#' @export
do.fisher.test <- function(fc.list, fg.ids, ml.list.up, ml.list.down, fc.threshold, ml.threshold=0.9, alternative = "two.sided"){
  
  if(length(fc.list) == 0){
    stop("fc.list is empty", immediate. = T)
  }
  
  fg.up_fc <- matrix(nrow=length(fc.list),ncol=dim(fc.list[[1]])[2],NA)
  fg.down_fc <- fg.up_fc
  
  bg.up_fc <- fg.up_fc
  bg.down_fc <- fg.up_fc
  
  bg.ids <- c()
  
  p.val.mat <- fg.up_fc
  for(l in 1:length(fc.list)){
    
    fg.ids <- fg.ids[fg.ids %in% rownames(ml.list.up[[l]])]
    fg.ids <- fg.ids[fg.ids %in% rownames(ml.list.down[[l]])]
    
    bg.ids <- rownames(ml.list.up[[l]])[! (rownames(ml.list.up[[l]]) %in% fg.ids)]
    bg.ids <- bg.ids[ bg.ids %in% rownames(ml.list.down[[l]])]
    
    fg.sub.up.ml <- na.omit(ml.list.up[[l]][fg.ids,]) > ml.threshold
    fg.sub.down.ml <- na.omit(ml.list.down[[l]][fg.ids,]) > ml.threshold
    
    bg.sub.up.ml <- na.omit(ml.list.up[[l]][bg.ids,]) > ml.threshold
    bg.sub.down.ml <- na.omit(ml.list.down[[l]][bg.ids,]) > ml.threshold
    
    fg.sub.fc.up.list <- fc.list[[l]][fg.ids,] > fc.threshold
    fg.sub.fc.down.list <- fc.list[[l]][fg.ids,] < -fc.threshold
    
    bg.sub.fc.up.list <- fc.list[[l]][bg.ids,] > fc.threshold
    bg.sub.fc.down.list <- fc.list[[l]][bg.ids,] < -fc.threshold
    
    fg.sub.up.mat <- colSums(fg.sub.up.ml & fg.sub.fc.up.list)
    fg.sub.down.mat <- colSums(fg.sub.down.ml & fg.sub.fc.down.list)
    
    bg.sub.up.mat <- colSums(bg.sub.up.ml & bg.sub.fc.up.list)
    bg.sub.down.mat <- colSums(bg.sub.down.ml & bg.sub.fc.down.list)
    
    fg.up_fc[l,] <- fg.sub.up.mat
    fg.down_fc[l,] <- fg.sub.down.mat
    
    bg.up_fc[l,] <- bg.sub.up.mat
    bg.down_fc[l,] <- bg.sub.down.mat
    
    rownames(fg.up_fc) <- rownames(bg.up_fc) <- names(ml.list.up)
    rownames(fg.down_fc) <- rownames(bg.down_fc) <- names(ml.list.down)
    
    for(j in 1:ncol(fg.up_fc)){
      N.table <- matrix(nrow=2,ncol=2,c(fg.up_fc[l,j], fg.down_fc[l,j], bg.up_fc[l,j], bg.down_fc[l,j]),byrow = F)
      p.val.mat[l,j] <- fisher.test(N.table,simulate.p.value = F, alternative = alternative)$p.value
    }
  }
  
  return(list(p.val.mat=p.val.mat, fg.up_fc=fg.up_fc, fg.down_fc=fg.down_fc, fg.ids = fg.ids))
}
#' @title Plot dendrogram colored by sample definition
#' @description Does a hierarchical clustering of the data and draws a dendrogam with colored branches and leaved based on the sample definition.
#' @author Alexander Gabel
#' @usage Compute_GO_Enrichment(geneUniverse, selectedGeneIds, pvalueCutoff=0.05, ontology="BP")
#' @param geneUniverse \code{vector} containing all gene ids of the whole dataset
#' @param selectedGeneIds \code{vector} containing subset of gene ids to check if their GO terms are enriched
#' @param pvalueCutoff numeric value defining the p-value cutoff for the adjusted p-values to call a GO term significantly enriched (default: 0.05)
#' @param ontology character string defining the GO ontology and must be one of "BP", "CC", or "MF"
#' @return a \code{list} 
#' @import org.At.tair.db
#' @import AnnotationDbi
#' @import GSEABase
#' @import Category
#' @export
Compute_GO_Enrichment <- function(geneUniverse, selectedGeneIds, pvalueCutoff=0.05, ontology="BP"){
  
  # The following code based on 
  # http://www.cbs.dtu.dk/chipcourse/Exercises/Ex_GO/GOexercise11.php

  geneUniverse <- toupper(geneUniverse)
  selectedGeneIds <- toupper(selectedGeneIds)
  
  frame <- toTable(org.At.tairGO)
  
  goframeData <- data.frame(frame$go_id, frame$Evidence, frame$gene_id)
  goFrame <- GOFrame(goframeData, organism = "Arabidopsis thaliana")
  
  suppressWarnings(goAllFrame <- GOAllFrame(goFrame))
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  params <- GSEAGOHyperGParams("Arabidopsis thaliana GO", geneSetCollection = gsc, geneIds = selectedGeneIds, universeGeneIds = geneUniverse, 
                                          ontology = ontology, pvalueCutoff = 1, conditional = FALSE, testDirection = "over")
  
  over <- GOstats::hyperGTest(params)
  over.df <- as.data.frame(summary(over))
  
  p_corrected <- p.adjust(over.df$Pvalue, method = "bonferroni")
  
  over.df <- data.frame(over.df[,1:2], Pvalue.adjusted = p_corrected, over.df[,-(1:2)])
  over.df <- over.df[order(over.df$'Pvalue.adjusted'),]
  
  return(over.df[over.df$'Pvalue.adjusted' < pvalueCutoff,])
  
}

#' @title Adjust list of p-value matrices
#' @description Adjusts all p-values in the list and separates them finally into the same matrix structure
#' @author Alexander Gabel
#' @usage adjust_and_split(p_val_list, adjust_method="BY")
#' @param p_val_list list of matrices containing uncorrected p-vales
#' @param adjust_method method for p-value correction. Default: "BY". For details see \link{p.adjust}
#' @return a \code{list} 
#' @export
adjust_and_split <- function(p_val_list, adjust_method="BY"){
  
  if(!is.list(p_val_list)){
    stop("Please check p_val_list. It has to be a list.")
  }
  
  if(!is.matrix(p_val_list[[1]])){
    stop("Please check p_val_list. Its elements have to be matrices.")
  }
  
  mat_rows <- nrow(p_val_list[[1]])
  mat_cols <- ncol(p_val_list[[1]])
  
  adjusted_p_vals_vec <- p.adjust(unlist(p_val_list), method = adjust_method)
  adjusted_p_vals_mat <- matrix(nrow=mat_rows, adjusted_p_vals_vec)
  
  j_mat <- 1
  adjusted_p_mat_list <- lapply(seq(1, ncol(adjusted_p_vals_mat)-mat_cols+1, mat_cols), 
                                function(i){ 
                                  mat <- adjusted_p_vals_mat[,i:(i+mat_cols-1)]
                                  rownames(mat) <- rownames(p_val_list[[j_mat]])
                                  colnames(mat) <- colnames(p_val_list[[j_mat]])
                                  j_mat <<- j_mat + 1
                                  mat
                                  }
                               )
  names(adjusted_p_mat_list) <- names(p_val_list)
  return(adjusted_p_mat_list)
}