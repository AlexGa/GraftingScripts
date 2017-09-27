#' @title Plot dendrogram colored by sample definition
#' @description Does a hierarchical clustering of the data and draws a dendrogam with colored branches and leaved based on the sample definition.
#' @author Alexander Gabel
#' @usage plot_colored_dendrogram(exp_data,groups, method="pearson", link="average", do.legend=F, unClusteredBranchColor="black", legend.pos = "topright", ...)
#' @param exp_data expression matrix, columns represent the samples and rows the genes (or transcripts)
#' @param groups numeric vector describing the relationships of columns to certain groups (like biological replicates)
#' @param method character string defining which correlation coefficient should be used as similarity measure pearson (default), spearman, or kendall
#' @param link character sting defining which agglomeration method to be used. Default: "average". For more information see parameter method at \link{hclust}
#' @param do.legend boolean value defining if a color legend should be ploted
#' @param unClusteredBranchColor character string defining the color to use if two brnaches belong to different sample types 
#' @param legend.pos character string defining the position of the legend e.g. "topright" (default), "topleft",... \link{legend}
#' @param ... optional plot parameters \link{plot}
#' @return a \code{plot} 
#' @export
plot_colored_dendrogram <- function(exp_data, groups, method="pearson", link="average", do.legend=F, unClusteredBranchColor="black", legend.pos = "topright", ...){
  
  if(!is.factor(groups)){
    groups <- factor(groups)
  }
  
  group_levels <- droplevels(groups)
  rep_indices <- split(seq_along(groups), f = group_levels)
  
  nCols <- length(unique(group_levels))
  
  labelColors <- c()
  
  if(nCols <= 10){
    labelColors <- suppressWarnings(RColorBrewer::brewer.pal(name = "Set1",n = nCols))
  }else{
    labelColors <- rainbow(nCols)
  }
  
  j <- 0
  repColors <- lapply(rep_indices, function(i){
    j <<- j +1
    cbind(colnames(exp_data)[i],labelColors[j])
  })
  
  sample.labels.colors <- do.call("rbind",repColors)
  rownames(sample.labels.colors) <- sample.labels.colors[,1]
  
  sample.labels <- sample.labels.colors[,1]
  
  dt.cor <- cor(exp_data,method = method)
  hc <- hclust(d = 1-as.dist(dt.cor),method = link)
  
  hcd = as.dendrogram(hc)
  
  # function to get color labels
  prevLeaf <- F
  prevNode <- c()
  test <- c()
  reachedFirstleaf <- F
  leaf.labels=NULL
  
  colLab <- function(n) {
    
    if (is.leaf(n)) {
      
      a <- attributes(n)
      labCol <- as.character(sample.labels.colors[a$label,2])
      if(!is.null(leaf.labels)){
        leafLab <- as.character(sample.labels.colors[a$label,1])
        attr(n, "label") <- leafLab
      }else{
        attr(n, "label") <- c(a$label)
      }
      attr(n, "edgePar") <- list(col = labCol,lwd=3);
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol,pch=NA,col=labCol)
      prevLeaf <<- T
      prevCol <<- labCol
      if(!reachedFirstleaf){
        reachedFirstleaf <<- T
      }
      prevNode <- n
    }
    else{
      if(reachedFirstleaf){
        
        if(attr(n[[1]],"edgePar")[1] == attr(prevNode,"edgePar")[1]){
          
        }
      }
      prevNode <- n
      if(reachedFirstleaf){
        dendrapply(n,)
      }
    }
    n
  }
  
  
  colTree <- function(n) {
    
    if (is.leaf(n)) {
      a <- attributes(n)
      #print(a)
      labCol <- as.character(sample.labels.colors[a$label,2])
      if(!is.null(leaf.labels)){
        leafLab <- as.character(sample.labels.colors[a$label,1])
        attr(n, "label") <- leafLab
      }else{
        attr(n, "label") <- c(a$label)
      }
      attr(n, "edgePar") <- list(col = labCol,lwd=3);
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol,pch=NA,col=labCol)
    }
    else{
      if(is.leaf(n[[1]]) & is.leaf(n[[2]])){
        # print("1")
        a <- attributes(n[[1]])
        labCol.1 <- as.character(sample.labels.colors[a$label,2])
        
        a <- attributes(n[[2]])
        labCol.2 <- as.character(sample.labels.colors[a$label,2])
        
        if(labCol.1 == labCol.2){
          attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
        }
      }
      if(is.leaf(n[[1]]) & !is.leaf(n[[2]])){
        
        a.1 <- attributes(n[[1]])
        # print(paste("2.1:  ",a.1$label))
        labCol.1 <- as.character(sample.labels.colors[a.1$label,2])
        
        a.2 <- attributes(n[[2]])
        labCol.2 <- a.2$edgePar$col
        
        # print("2.2")
        attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
        
        if(!is.null(labCol.2)){
          
          if(labCol.1 == labCol.2){
            attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
          }
        }
      }
      if(!is.leaf(n[[1]]) & is.leaf(n[[2]])){ ## Fall tritt niemals ein, da Baum LWR-Traversiert wird
        # print("3")
        a.2 <- attributes(n[[2]])
        labCol.2 <- as.character(sample.labels.colors[a.2$label,2])
        
        a.1 <- attributes(n[[1]])
        labCol.1 <- a.1$edgePar$col
        
        #attr(n, "edgePar") <- list(col = labCol.2,lwd=3);
        #print(paste0(labCol.1," -->",labCol.2))
        if(!is.null(labCol.1)){
          
          if(labCol.1 == labCol.2){
            
            attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
          }
        }
      }
      if(!is.leaf(n[[1]]) & !is.leaf(n[[2]])){
        # print("4")
        a.1 <- attributes(n[[1]])
        labCol.1 <- a.1$edgePar$col
        
        # print(labCol.1)
        
        if(!is.null(labCol.1)){
          a.2 <- attributes(n[[2]])
          labCol.2 <- a.1$edgePar$col
          if(!is.null(labCol.2)){
            if(labCol.1 == labCol.2){
              attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
              attr(n[[2]], "edgePar") <- list(col = labCol.1,lwd=3);
            }
          }
        }
      }
    }
    n
  }
  # using dendrapply
  
  clusDendro = dendrapply(hcd, colTree)
  
  rownames(sample.labels.colors) <- sample.labels.colors[,1]
  names(sample.labels) <- sample.labels.colors[,1]
  clusDendro = dendrapply(clusDendro, colTree)
  
  clean.dendro <- function(n){
    
    if(!is.leaf(n)){
      a.1 <- attributes(n[[1]])
      # print(paste("2.1:  ",a.1$label))
      labCol.1 <- a.1$edgePar$col
      
      a.2 <- attributes(n[[2]])
      labCol.2 <- a.2$edgePar$col
      if(!is.null(labCol.1) & !is.null(labCol.2)){
        if(labCol.1 != labCol.2){
          attr(n, "edgePar") <- list(col = unClusteredBranchColor,lwd=3);
        }
        if(labCol.1 == labCol.2){
          attr(n, "edgePar") <- list(col = labCol.1,lwd=3);
        }
      }
      
    }
    return(n)
  }
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  clusDendro = dendrapply(clusDendro, clean.dendro)
  
  plot(as.dendrogram(clusDendro,hang = -1,check = F),...)
  if(do.legend){
    legend(legend.pos,legend = legend.text,bty="n",fill=labelColors,border = NA,ncol = 2)
  }
  
}


#' @title Principle component analysis (PCA) or multidimensional scaling plot between samples
#' @description Plot samples on a two-dimensional scatterplot in which the distances on the plot indicate the expression differences between the samples
#' @author Alexander Gabel
#' @usage plotPCA(exp_data, groups, do.legend=F, plot_time_points=F, log=T, do.MDS=F, cols=NULL, epsilon = 1, return_plotMatrix = FALSE, ...)
#' @param exp_data expression matrix, columns represent the samples and rows the genes (or transcripts)
#' @param groups numeric vector describing the relationships of columns to certain groups (like biological replicates)
#' @param do.legend boolean value defining if a color legend should be ploted
#' @param plot_time_points boolean value indicating if the grafting time points shall be plotted
#' @param log boolean value defining if the expression data should be log transformed before the calculation of PCA or MDS
#' @param do.MDS boolean value defines if PCA or MDS should be calculated. Default: TRUE
#' @param cols character vector defining the colors to use for the different samples
#' @param epsilon numeric value indicates pseudocount used if \code{log=TRUE} 
#' @param return_plotMatrix boolean value indicates if PCA or MDS results should be returned 
#' @param ... optional plot parameters \link{plot}
#' @return a \code{matrix} 
#' @export
plotPCA <- function(exp_data, groups, do.legend=F, plot_time_points=F, log=T, do.MDS=F, cols=NULL, epsilon = 1, return_plotMatrix = FALSE, ...){
  
  if(!is.matrix(exp_data) && !is.data.frame(exp_data)){
    warning("exp_data has to be a matrix or data.frame.")
    return()
  }
  
  nSamples <- length(unique(groups))
  # Expression matrix exp_data -> samples represented by columns
  exp_data <- t(as.matrix(exp_data))
  
  if(log){
    exp_data <- log(exp_data + epsilon)
  }
  
  scaling_mat <- c()
  x <- y <- c()
  
  if(do.MDS){
    
    exp_dist <- dist(exp_data)
    scaling_mat <- cmdscale(exp_dist, eig=TRUE, k=2)
    x <- scaling_mat$points[,1]
    y <- scaling_mat$points[,2]
  }else{
    
    scaling_mat <- prcomp(exp_data)
    x <- scaling_mat$x[,1]
    y <- scaling_mat$x[,2]
  }
  
  if(is.null(cols)){
    if(nSamples > 9){
      cols <- rainbow(n = nSamples, alpha = 0.8)
    }else{
      cols <- RColorBrewer::brewer.pal(name = "Set1",n = nSamples)
      cols <- scales::alpha(cols,alpha = 0.8) 
    }
  }
  
  gr_table <- table(groups)
  colNumbs <- unlist(sapply(1:nSamples,function(i){
    rep(i,each=gr_table[i])
  }))
  
  time_points <- unlist(lapply(strsplit(x = rownames(exp_data),split = " "),function(i)(i[1])))
  
  if(do.MDS){
    plot(x, y, col = cols[colNumbs], xlab = "Coordinate 1", ylab = "Coordinate 2", pch=19,...)
  }else{
    plot(x, y, col = cols[colNumbs], xlab = paste0("PC1 (", round(scaling_mat$sdev[1], digits = 2),"%)"), ylab = paste0("PC2 (",round(scaling_mat$sdev[2], digits = 2),"%)"), pch=19,...)
  }
  
  if(plot_time_points){
    text(x = x, y=y, labels = time_points)
  }
  
  legend.text <- c()
  
  if(do.legend){
    legend.text <- unique(unlist(lapply(strsplit(rownames(exp_data), split = " "), function(i)paste0(i[-1], collapse=" "))))
    legend("topleft", legend = legend.text, col = cols[1:nSamples], bty = "n", pch = 19, ncol=2, pt.cex = 2)    
  }
  
  if(return_plotMatrix){
    return(scaling_mat)
  }
}

#' @title Barplot of up- and down- regulated differentially expressed genes 
#' @description Plots the relative or absolute number of up- and down- regulated genes and a asterisks if the ratio is significant 
#' @author Alexander Gabel
#' @usage barplot_up_down(up_mat, down_mat, p_val_mat, names.arg, main, labels)
#' @param up_mat,down_mat (relative) number of up or down regulated genes
#' @param p_val_mat matrix of p-values for each time point and each treatment. Indicating if ratio of up- and down- regulated genes is significant compared to the background.
#' @param ylim defines the limits of the y-axis
#' @param names.arg defines the names of the bars. Details: \link{barplot}
#' @param main defines the title of the plot
#' @param labels defines the labels of the x-axis
#' @param ylab defines the label of the y-axis
#' @return a barplot 
#' @export
barplot_up_down <- function(up_mat, down_mat, p_val_mat, ylim = c(-1.2,1.45), names.arg, main, labels, ylab="Proportion of overlap"){
  
  bar.coords <- barplot(up_mat, col = "green4", beside = T, las = 2, ylim = ylim, main = main, names.arg = names.arg, ylab = "")
  barplot(down_mat, col="red4", beside=T, add=T, yaxt="n", xaxt="n")
  
  mtext(text = ylab, side = 2, line = 2.75, at = 0, cex=1.5)
  text(x = bar.coords[4,], y=1.28, labels = labels)
  
  star.line.coords <- apply(up_mat, 2, function(up_row){
    line_coords <- c()
    if(max(up_row) + 0.1 > 0.5){
      line_coords <- rep(max(up_row) + 0.1, length(up_row))
    }else{
      line_coords <- rep(0.5, length(up_row))
    }
  })
  
  if(sum(t(p_val_mat) < 5e-2) > 0){
    text(x = bar.coords[t(p_val_mat) < 5e-2], y = star.line.coords[t(p_val_mat) < 5e-2], labels = "*")
  }
  
  abline(h = 0)
  vert_line_coords <- as.vector(bar.coords[c(1,8),])
  vert_line_coords <- vert_line_coords[c(-1,-length(vert_line_coords))]
  vert_line_coords <- matrix(data = vert_line_coords,nrow = 2,ncol = 3,byrow = F)
  vert_line_coords <- apply(vert_line_coords,2,function(i){
    i[1] + (diff(i)/2)
  })
  abline(v = vert_line_coords,lty=2)
  legend("bottomleft",legend = c("DEG up", "DEG down"),fill=c("green4","red4"),cex = 0.75)
}

#' @title Barplot of up- and down- regulated differentially expressed genes 
#' @description Plots the relative or absolute number of up- and down- regulated genes and a asterisks if the ratio is significant 
#' @author Alexander Gabel
#' @usage barplot_graft_formation(count_mat, pval_mat, ylab="# genes intersecting", pval_threshold = 0.05, main)
#' @param count_mat number of genes during each time point for each condition
#' @param pval_mat matrix of p-values for each time point and each treatment. Indicating if ratio of up- and down- regulated genes is significant compared to the background.
#' @param pval_threshold p-value threshold to indocate if a p-value is significant 
#' @param main defines the title of the plot
#' @param ylab defines the label of the y-axis
#' @return a barplot 
#' @export
barplot_graft_formation <- function(count_mat, pval_mat, ylab="# genes intersecting", pval_threshold = 0.05, main){
  
  ncols <- RColorBrewer::brewer.pal(n = nrow(count_mat), name="Set1")[c(3,1,2)]
  y_lim <- c(0, max(count_mat)+max(count_mat)*0.25)
  
  b_plt <- barplot(count_mat, ylim=y_lim, beside = T, col = ncols, names = gsub(pattern = "_.+", replacement = "", x = colnames(count_mat)), ylab = ylab, main = main)
  
  legend("topleft",legend = gsub(rownames(count_mat), pattern = "_", replacement = " & "), bty="n", fill = ncols)
  
  
  if(sum(pval_mat < pval_threshold) > 0){
    text(x = b_plt[pval_mat < pval_threshold], y=count_mat[pval_mat < pval_threshold], labels = "*", pos = 3) 
  }
  
}
