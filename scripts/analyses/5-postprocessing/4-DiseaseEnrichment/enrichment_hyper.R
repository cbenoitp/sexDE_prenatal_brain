#------------------------------------------------------------------#
# Functions to do the enrichment analysis with hypergeometric test #
#------------------------------------------------------------------#


#### Function to compute the enrichment and depletion p-value with the hypergeometric test ####
enrichment_stat <- function(overlap_nbr, geneList_nbr, DE_nbr, universe_nbr){
  # overall enrichment
  #print(">> Enrichment")
  pvalEnrich <- phyper(overlap_nbr-1, geneList_nbr, universe_nbr - geneList_nbr, DE_nbr, lower.tail = FALSE)
  # overall depletion
  #print(">> Depletion")
  pvalDeple <-phyper(overlap_nbr,  geneList_nbr, universe_nbr -  geneList_nbr, DE_nbr, lower.tail = TRUE)
  # effect size
  #print(">> Relative enrichment")
  relatEnrich <- (overlap_nbr/ geneList_nbr)/(DE_nbr/universe_nbr)
  
  return(list("pvalEnrich" = pvalEnrich, "pvalDeple" = pvalDeple, "relatEnrich" = relatEnrich))
}



#### Function to select the signif DE genes in the different analysis ####
select_DEsignif <- function(DEresults, pval_threshold = 0.05, direction = "both", bayes = FALSE, bayesDataType = NULL){
  if (bayes){
    DEresults_signif <- DEresults[which(DEresults$specificity == bayesDataType), ]
    if (direction == "up"){
      DEresults_signif <- DEresults_signif[which(DEresults_signif$logFC_prenatal > 0), ]
    } else if (direction == "down"){
      DEresults_signif <- DEresults_signif[which(DEresults_signif$logFC_prenatal < 0), ]
    }
  } else{
    if (direction == "both"){
      DEresults_signif <- DEresults[which(DEresults$adj.P.Val < pval_threshold), ]
    } else if (direction == "up"){
      DEresults_signif <- DEresults[which(DEresults$adj.P.Val < pval_threshold & DEresults$logFC > 0), ]
    } else if (direction == "down"){
      DEresults_signif <- DEresults[which(DEresults$adj.P.Val < pval_threshold & DEresults$logFC < 0), ]
    }
  }
  
  return(DEresults_signif)
}



#### Function to select the random gene sets for the "permutation" analysis ####
generateRandomDE <- function(DE, DEsignif, ctrlExpr = FALSE){
  if (ctrlExpr){
    nbrDE_byExprBin <- as.vector(table(DEsignif$exprbin))
    # select a random set of genes for each bin of expression
    randomDEgeneset_idx <- c()
    for (binExpr in 1:10){
      bin_genes_idx <- which(DE$exprbin == binExpr)
      randomDEgeneset_idx_bin <- sample(bin_genes_idx, nbrDE_byExprBin[binExpr])
      randomDEgeneset_idx <- c(randomDEgeneset_idx, randomDEgeneset_idx_bin)
    }
  } else{
    randomDEgeneset_idx <- sample(1:length(DE$gene_name), length(DEsignif$gene_name)) 
  }
  randomDE <- DE[randomDEgeneset_idx, ]
  return(randomDE)
}



#### Function to compute permutation p-value ####
compute_permPvalue <- function(value, randomValues){
  nbrUp <- sum(randomValues > value)
  nbrLow <- sum(randomValues < value)
  if (value > 1){
    pval <- nbrUp / length(randomValues)
  } else if (value < 1){
    pval <- nbrLow / length(randomValues)
  }
  
  return(pval)
}



#### Function to compute the enrichment/depletion of a specific geneList in DE analysis results ####
enrichment_analysis <- function(DEresults, DEresults_signif, geneList, analysisName, geneListName, pval_threshold = 0.05, 
                                direction = "both", bayes = FALSE, bayesDataType = NULL, geneName = TRUE){
  universe_nbr <- length(DEresults$gene_name)
  DE_nbr <- length(DEresults_signif$gene_name)
  # for the geneList, keep only gene that are in the DEresults
  if (geneName){
    geneList_filtered <- geneList[which(geneList %in% DEresults$gene_name)]
  } else{
    geneList_filtered <- geneList[which(geneList %in% DEresults$ensemblID)]
  }
  geneList_nbr <- length(geneList_filtered)
  # overlap
  if (geneName){
    overlap_nbr <- length(intersect(DEresults_signif$gene_name, geneList_filtered))
  } else{
    overlap_nbr <- length(intersect(DEresults_signif$ensemblID, geneList_filtered))
  }
  # get enrichment stat and effect size
  enrichment_stat_res <- enrichment_stat(overlap_nbr, geneList_nbr, DE_nbr, universe_nbr)
  # format data
  results <- data.frame("analysis" = analysisName, "direction" = direction, "geneList" = geneListName,
                        "nbr_DE" = DE_nbr, "nbr_geneList" = geneList_nbr, 
                        "overlap" = overlap_nbr, "N" = universe_nbr, 
                        "pval_enrich" = enrichment_stat_res$pvalEnrich, "pval_deple" = enrichment_stat_res$pvalDeple,
                        "relative_enrich" = enrichment_stat_res$relatEnrich)
  return(results)
}



#### Function to compute the enrichment/depletion of a list of gene sets in DE analysis results ####
#### + compare this enrichment/depletion to a number of random gene sets                        ####
enrichment_analysis_permut <- function(DEresults, geneList, analysisName, geneListNames, nbr_permut = 1000, pval_threshold = 0.05, 
                                       direction = "both", bayes = FALSE, bayesDataType = NULL, ctrlExpr = FALSE, geneName = TRUE){
  # select signif DE 
  DEresults_signif <- select_DEsignif(DEresults = DEresults, pval_threshold = pval_threshold, direction = direction, bayes = bayes, bayesDataType = bayesDataType)
  # for each gene sets in the geneList
  results_enrich_geneset_full <- data.frame("analysis" = NA, "direction" = NA, "geneList" = NA, "nbr_DE" = NA, "nbr_geneList" = NA, 
                                            "overlap" = NA, "N" = NA, "pval_enrich" = NA, "pval_deple" = NA, "relative_enrich" = NA)
  results_enrich_geneset_full <- results_enrich_geneset_full[0,]
  for (i in 1:length(geneList)){
    geneset <- geneList[[i]]
    geneListName <- geneListNames[[i]]
    # enrichment for geneset
    results_enrich_geneset <- enrichment_analysis(DEresults = DEresults, DEresults_signif = DEresults_signif, geneList = geneset, 
                                                  geneListName = geneListName, analysisName = analysisName, pval_threshold = pval_threshold,
                                                  direction = direction, bayes = bayes, bayesDataType = bayesDataType, geneName = geneName)
    results_enrich_geneset_full <- rbind(results_enrich_geneset_full, results_enrich_geneset)
  }
  # enrichment for 'nbr_permut" random genesets of the same size of the DE signif
  results_enrich_randomGeneset_full <- results_enrich_geneset_full[0,]
  for (i in 1:nbr_permut){
    # get a random geneset of the same size as the DEresults_signif
    random_DE <- generateRandomDE(DEresults, DEresults_signif, ctrlExpr = ctrlExpr)
    # analysis done for every geneList
    for (j in 1:length(geneList)){
      geneset <- geneList[[j]]
      geneListName <- geneListNames[[j]]
      # enrichment for random geneset
      results_enrich_randomGeneset <- enrichment_analysis(DEresults = DEresults, DEresults_signif = random_DE, geneList = geneset, 
                                                          geneListName = geneListName, analysisName = analysisName, pval_threshold = pval_threshold,
                                                          direction = direction, bayes = bayes, bayesDataType = bayesDataType, geneName = geneName)
      results_enrich_randomGeneset_full <- rbind(results_enrich_randomGeneset_full, results_enrich_randomGeneset)
    }
  }
  # for each gene set in geneList, merge the results from the random genesets and add them to final results
  results_enrich_geneset_full$pval_enrich_random <- NA
  results_enrich_geneset_full$pval_deple_random <- NA
  results_enrich_geneset_full$relative_enrich_random <- NA
  results_enrich_geneset_full$permut_pval <- NA
  for (geneset in geneListNames){
    tmp_results_enrich <- results_enrich_randomGeneset_full[which(results_enrich_randomGeneset_full$geneList == geneset),]
    # add info to final results
    results_enrich_geneset_full$pval_enrich_random[which(results_enrich_geneset_full$geneList == geneset)] <- paste(tmp_results_enrich$pval_enrich, sep = ";", collapse = ";")
    results_enrich_geneset_full$pval_deple_random[which(results_enrich_geneset_full$geneList == geneset)] <- paste(tmp_results_enrich$pval_deple, sep = ";", collapse = ";")
    results_enrich_geneset_full$relative_enrich_random[which(results_enrich_geneset_full$geneList == geneset)] <- paste(tmp_results_enrich$relative_enrich, sep = ";", collapse = ";")
    results_enrich_geneset_full$permut_pval[which(results_enrich_geneset_full$geneList == geneset)] <- compute_permPvalue(results_enrich_geneset_full$relative_enrich[which(results_enrich_geneset_full$geneList == geneset)], tmp_results_enrich$relative_enrich)
  }
  results_enrich_geneset_full$ctrlExpr <- ctrlExpr
  
  return(results_enrich_geneset_full)
}


