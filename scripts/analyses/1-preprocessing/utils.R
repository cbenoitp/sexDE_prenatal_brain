# functions to help analysis

#' Function to load a gtf file.
#' @param GTFfile String containing the path and the name of the gtf file to load.
#' @param type String defining if the information of the transcripts or the genes should be loaded.
#' @return Data frame containing the annotation information.
#' @examples
#' 
loadGTFfile <- function(GTFfile, type = "transcript"){
  gtf_annot <- rtracklayer::import(GTFfile)
  annot <- gtf_annot[which(elementMetadata(gtf_annot)$type == type),]
  if (type == "transcript"){
    annot <- data.frame(transcript_id = elementMetadata(annot)$transcript_id,
                        gene_id = elementMetadata(annot)$gene_id,
                        transcript_name = elementMetadata(annot)$transcript_name,
                        gene_name = elementMetadata(annot)$gene_name,
                        transcript_type = elementMetadata(annot)$transcript_type,
                        gene_type = elementMetadata(annot)$gene_type,
                        chr = seqnames(annot))
  } else if (type == "gene"){
    annot <- data.frame(gene_id = elementMetadata(annot)$gene_id,
                        gene_name = elementMetadata(annot)$gene_name,
                        gene_type = elementMetadata(annot)$gene_type,
                        chr = seqnames(annot))
  }
  return(annot)
}


#' Function to assign sex depending on the expression of XIST and the Y-chr protein coding genes.
#' @param xist Vector of expression of the XIST gene.
#' @param chry Vector of the sum of expression of the Y-chr protein coding genes.
#' @param threshold Integer value defining the threshold to assign sex.
#' @return Vector of assigned sex (0 = Unassigned, 1 = Male & 2 = Female).
#' @examples
#' 
assignSexFromExpression <- function(xist, chry, threshold){
  assigned_sex <- ifelse(chry > threshold & xist < threshold, 
                         1, ifelse(chry_sum < threshold & xist > threshold, 2, 0))
  return(assigned_sex)
}


#' Function to plot the sex of the samples on the expression of XIST and Y-chr protein coding genes.
#' @param xist Vector of expression of the XIST gene.
#' @param chry Vector of the sum of expression of the Y-chr protein coding genes.
#' @param sex  Vector of sex (0 = Unassigned, 1 = Male & 2 = Female).
#' @param threshold Interger value defining the threshold for the assignment of sex.
#' @return Ggplot object containing the plot.
#' @examples
#' 
plotSexOnExpression <- function(xist, chry, sex, threshold){
  p <- ggplot()
  p <- p + geom_hline(yintercept = threshold)
  p <- p + geom_vline(xintercept = threshold)
  p <- p + geom_point(aes(x=xist[sex==0], y=chry[sex==0]), 
                      colour="#C9C9C9")
  p <- p + geom_point(aes(x=xist[sex==2], y=chry[sex==2]), 
                      colour="#c42134")
  p <- p + geom_point(aes(x=xist[sex==1], y=chry[sex==1]), 
                      colour="#074c9b")
  p <- p + xlab("Expression of XIST") + ylab("Expression of protein-coding genes on Y chr")
  p <- p + scale_x_log10() + scale_y_log10()
  p <- p + theme_bw()
  return(p)
}


#' Plot the individuals of a PCA in a plan
#' @param PCAresults Output object from ade4 dudi.pca function.
#' @param colorFactor Vector of factors. It should have the same length as the number of individuals.
#' @param factorName String corresponding to the name of factor used for coloring. It will be used as a legend title.
#' @param axis1 Integer indicating which principal component to use as x axis.
#' @param axis2 Integer indicating which principal component to use as y axis.
#' @param color Vector of color. Should be the same length as the number of levels in \code{colorFactor}. Default: NULL
#' @return Ggplot2 object containing the plot colored by the \code{colorFactor}.
PCAplot <- function(PCAresults, colorFactor, factorName, axis1 = 1, axis2 = 2, color = NULL) {
  # Compute the percent of variance explained by each eighen values
  varExplained <- PCAresults$eig/sum(PCAresults$eig)*100
  # Define the color if it is not already defined
  if (is.null(color)){
    color <- rainbow(length(levels(colorFactor)))
  }
  # Plot
  p <- fviz_pca_ind(PCAresults,
                    axes = c(axis1, axis2),
                    geom = "point",
                    habillage = colorFactor, # colorer par groupes
                    palette = color,
                    legend.title = factorName,
                    invisible = "quali")
  p <- p + labs(title = "PCA" , 
                x = sprintf("PC%s (%s %%)", axis1, round(varExplained[axis1])), 
                y = sprintf("PC%s (%s %%)", axis2, round(varExplained[axis2])))
  p <- p + theme_minimal()
  return(p)
}

#' Plot the individuals and the first few variables of a PCA in a plan
#' @param PCAresults Output object from ade4 dudi.pca function.
#' @param colorFactor Vector of factors. It should have the same length as the number of individuals.
#' @param factorName String corresponding to the name of factor used for coloring. It will be used as a legend title.
#' @param axis1 Integer indicating which principal component to use as x axis.
#' @param axis2 Integerindicating which principal component to use as y axis.
#' @param color Vector of color. Should be the same length as the number of levels in \code{colorFactor}. Default: NULL
#' @param nbrVar Integer indicating the number of the top variables to plot
#' @return Ggplot2 object containing the plot colored by the \code{colorFactor}.
PCAbiplot <- function(PCAresults, colorFactor, factorName, axis1 = 1, axis2 = 2, color = NULL, nbrVar = 10) {
  # Compute the percent of variance explained by each eighen values
  varExplained <- PCAresults$eig/sum(PCAresults$eig)*100
  # Define the color if it is not already defined
  if (is.null(color)){
    color <- rainbow(length(levels(colorFactor)))
  }
  # Plot
  p <- fviz_pca_biplot(PCAresults,
                       axes = c(axis1, axis2),
                       geom.ind = "point",
                       geon.var = "arrow",
                       habillage = colorFactor, # colorer par groupes
                       palette = color,
                       legend.title = factorName,
                       select.var = list(contrib = nbrVar),
                       col.var = "grey",
                       invisible = "quali")
  p <- p + labs(title = "PCA" , 
                x = sprintf("PC%s (%s %%)", axis1, round(varExplained[axis1])), 
                y = sprintf("PC%s (%s %%)", axis2, round(varExplained[axis2])))
  p <- p + theme_minimal()
  return(p)
}


#' Plot the number of male and female samples for a given covariate
#' @param sex Vector of 1 and 2 corresponding to the sex of the samples.
#' @param covariate A vector of factors. It should have the same length and the same order as the sex vector.
#' @param covariateName String corresponding to the name of factor of interest. It will be used as the x axis title.
#' @param covariate2 A vector of factors corresponding to a 2nd covariate. It will be used for the faceting of the plots.
#' It should have the same length and the same order as the sex vector.
#' @param covariate2Name String corresponding to the name of 2nd covariate of interest.
#' @return Return a ggplot2 object containing the plot of the repartition of sex for the covariate.
PlotSexRepartition <- function(sex, covariate, covariateName, covariate2 = NULL, covariate2Name = NULL, ylim = NULL, unknownSex = FALSE){
  # compute the number of male and female samples for each level of the covariate
  if (is.null(covariate2)){
    nbr <- table(covariate, sex)
    # reformat the data for the plot
    nbr_plot <- reshape2::melt(nbr)
    names(nbr_plot) <- c(covariateName, "sex", "nbr")
  } else{
    nbr_plot <- data.frame(covariate1=factor(),
                           covariate2=factor(), 
                           sex=factor(),
                           nbr=integer())
    for (cov in levels(covariate)){
      # add the number of samples for each sex in the different levels 
      nbr <- table(covariate2[covariate == cov], sex[covariate == cov])
      new_rows <- reshape2::melt(nbr)
      names(new_rows) <- c("covariate2", "sex", "nbr")
      new_rows$covariate <- as.factor(cov)
      nbr_plot <- rbind(nbr_plot, new_rows[, c(4,1,2,3)])
    }
    names(nbr_plot) <- c(covariateName, covariate2Name, "sex", "nbr")
  }
  # plot
  p <- ggplot(nbr_plot, aes(x=nbr_plot[,covariateName], y=nbr, group=as.factor(sex), color=as.factor(sex)))
  p <- p + geom_point() + geom_line() 
  if (!is.null(covariate2)){
    p <- p + facet_grid(nbr_plot[,covariate2Name]~.)
  }
  if (unknownSex){
    p <- p + scale_color_manual(values = c("0"="#C9C9C9", "1"="#074c9b", "2"="#c42134"), name = "Sex", breaks = c("0", "1", "2"), labels = c("Unknown", "Male", "Female"))
  }else{
    p <- p + scale_color_manual(values = c("1"="#074c9b", "2"="#c42134"), name = "Sex", breaks = c("1", "2"), labels = c("Male", "Female"))
  }
  if (!is.null(ylim)){
    p <- p + ylim(ylim)
  }
  p <- p + xlab(covariateName) + ylab("Number of samples") + theme_bw()
  p <- p + theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
  return(p)
}
