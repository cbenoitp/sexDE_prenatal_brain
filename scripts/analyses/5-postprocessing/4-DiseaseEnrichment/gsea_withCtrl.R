#-----------------------------------------#
# Functions to do the enrichment analysis #
#-----------------------------------------#

## Function gsea_withCtrl

gsea_withCtrl <- function(pathways, stats, nperm, p = 1, withCtrl = FALSE, ctrlVariable = NULL, 
                          scoreType = c("std", "pos", "neg"), printDebug = FALSE){
  if (printDebug){
    print("----> start gsea analysis")  
  }

  
  # 1- Prepare pathways and stats
  if (printDebug){
    print("--------> prepare data gsea")
  }
  scoreType <- match.arg(scoreType)
  pp <- prepare_gsea_data(pathways, stats, p = p, scoreType = scoreType, withCtrl = withCtrl, 
                          ctrlVariable = ctrlVariable, printDebug = printDebug)
  pathwaysFiltered <- pp$filtered
  pathwaysSizes <- pp$sizes
  stats <- pp$stats
  ctrlDistrib <- pp$ctrlDistrib
  
  # 3- For each pathway, compute ES
  if (printDebug){
    print("--------> compute ES")
  }
  gseaStatRes <- do.call(rbind, lapply(pathwaysFiltered, compute_ES,
                                       stats=stats,
                                       scoreType=scoreType,
                                       printDebug=printDebug))
  leadingEdges <- mapply("[", list(names(stats)), gseaStatRes[, "leadingEdge"], SIMPLIFY = FALSE)
  pathwayScores <- unlist(gseaStatRes[, "res"])
  
  # 4- Do the GSEA analysis
  # do a permutation test to estimate the nominal p-value,
  # compute the normalized ES (NES)
  # and apply the mutiple testing correction
  if (printDebug){
    print("--------> compute pvalues")
  }
  pvals <- compute_pval(pathwayScores, pathwaysSizes, pathwaysFiltered, 
                        leadingEdges, nperm, stats, scoreType,
                        withCtrl, ctrlVariable, ctrlDistrib, printDebug = printDebug)
  
  pvals[, nLeZero := NULL]
  pvals[, nGeZero := NULL]
  pvals[, leZeroMean := NULL]
  pvals[, geZeroMean := NULL]
  pvals[, nLeEs := NULL]
  pvals[, nGeEs := NULL]
  
  setcolorder(pvals, c("pathway", "pval", "padj", "ES", "NES",
                       "nMoreExtreme", "size", "leadingEdge"))
  # Makes pvals object printable immediatly
  pvals <- pvals[]
  
  return(pvals)
}


## Function prepare_gsea_data

prepare_gsea_data <- function(pathways, stats, p = 1, withCtrl = FALSE, ctrlVariable = NULL, scoreType = "std", 
                              printDebug = FALSE){
  # Sort and prepare the stats vector
  if (scoreType %in% c("pos", "neg")){
    stats <- abs(stats)
  }
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  #stats <- sort(stats, decreasing = TRUE)
  stats <- stats[ord]
  stats <- sign(stats) * (abs(stats) ^ p)
  stats <- stats / max(abs(stats))
  
  # Prepare the pathways
  pathwaysFiltered <- lapply(pathways, function(x) {sort(unique(na.omit(fmatch(x, names(stats))))) })
  pathwaysSizes <- sapply(pathwaysFiltered, length)
  
  toKeep <- which(pathwaysSizes > 0)
  pathwaysFiltered <- pathwaysFiltered[toKeep]
  pathwaysSizes <- pathwaysSizes[toKeep]
  
  # Result
  res <- list(filtered = pathwaysFiltered,
              sizes = pathwaysSizes,
              stats = stats)
  
  # Distributions of the control variable
  if (withCtrl){
    ctrlDistrib <- lapply(pathwaysFiltered, function(x) {ctrlVariable[x] })
    ctrlDistrib <- lapply(ctrlDistrib, table)
    res[["ctrlDistrib"]] <- ctrlDistrib
  }
  
  return(res)
}


## Function compute_ES

compute_ES <- function(stats, selectedStats, p = 1, returnAllExtremes = FALSE, scoreType = c("std", "pos", "neg"),
                       printDebug = FALSE){
  scoreType <- match.arg(scoreType)
  # this function is mostly a copy of the function "calcGseaStat" from the fgsea package
  S <- selectedStats
  r <- stats
  
  S <- sort(S)
  
  m <- length(S)
  N <- length(r)
  NR <- (sum(abs(r[S])^p))
  rAdj <- abs(r[S]) ^ p
  
  if (NR == 0) {
    rCumSum <- seq_along(rAdj) / length(rAdj)
  } else {
    rCumSum <- cumsum(rAdj) / NR
  }
  
  tops <- rCumSum - (S - seq_along(S)) / (N - m)
  if (NR == 0) {
    bottoms <- tops - 1 / m
  } else {
    bottoms <- tops - rAdj / NR
  }
  
  if (printDebug){
    browser()
  }
  
  maxP <- max(tops)
  minP <- min(bottoms)
  
  switch(scoreType,
         std = geneSetStatistic <- ifelse(maxP == -minP, 0, ifelse(maxP > -minP, maxP, minP)),
         pos = geneSetStatistic <- maxP,
         neg = geneSetStatistic <- minP)
  
  res <- list(res=geneSetStatistic)
  if (returnAllExtremes) {
    res <- c(res, list(tops=tops, bottoms=bottoms))
  }
  leadingEdge <- if (maxP > -minP) {
    S[seq_along(S) <= which.max(bottoms)]
  } else if (maxP < -minP) {
    rev(S[seq_along(S) >= which.min(bottoms)])
  } else {
    NULL
  }
  res <- c(res, list(leadingEdge=leadingEdge))
  
  return(res)
}


## Function compute_pval

compute_pval <- function(pathwayScores, pathwaysSizes, pathwaysFiltered, leadingEdges, nperm, stats, scoreType,
                         withCtrl = FALSE, ctrlVariable = NULL, ctrlDistrib = NULL, printDebug = FALSE){
  
  universe <- seq_along(stats)
  nbrPathways <- length(pathwaysSizes)
  counts <- list()
  
  # for each permutation
  for (i in seq_len(nperm)) {
    
    leEs <- rep(0, nbrPathways)
    geEs <- rep(0, nbrPathways)
    leZero <- rep(0, nbrPathways)
    geZero <- rep(0, nbrPathways)
    leZeroSum <- rep(0, nbrPathways)
    geZeroSum <- rep(0, nbrPathways)
    
    randEsP <- rep(0, nbrPathways)
    
    # for each pathway
    for (k in seq_len(nbrPathways)){
      
      if (printDebug){
        browser()
      }
      
      randSample <- sort(randomSelection(length(universe), pathwaysSizes[k], withCtrl, ctrlVariable, ctrlDistrib[[k]],
                                         printDebug = printDebug))
      randEsP_k <- compute_ES(
        stats = stats,
        selectedStats = randSample,
        scoreType = scoreType,
        p = 1,
        printDebug = printDebug)
      randEsP[k] <- randEsP_k$res
    }
    leEs <- leEs + (randEsP <= pathwayScores)
    geEs <- geEs + (randEsP >= pathwayScores)
    leZero <- leZero + (randEsP <= 0)
    geZero <- geZero + (randEsP >= 0)
    leZeroSum <- leZeroSum + pmin(randEsP, 0)
    geZeroSum <- geZeroSum + pmax(randEsP, 0)
    
    counts[[i]] <- data.table(pathway=seq_len(nbrPathways),
                              leEs=leEs, geEs=geEs,
                              leZero=leZero, geZero=geZero,
                              leZeroSum=leZeroSum, geZeroSum=geZeroSum)
  }
  counts <- rbindlist(counts)
  
  pvals <- counts[, list(leZeroMean = sum(leZeroSum, na.rm = T) / sum(leZero, na.rm = T),
                         geZeroMean = sum(geZeroSum, na.rm = T) / sum(geZero, na.rm = T),
                         nLeZero = sum(leZero, na.rm = T),
                         nGeZero = sum(geZero, na.rm = T),
                         nLeEs = sum(leEs, na.rm = T),
                         nGeEs = sum(geEs, na.rm = T)),
                  by = .(pathway)]
  
  pvals[, ES := pathwayScores[pathway]]
  pvals[, NES := as.numeric(NA)]
  
  switch(scoreType,
         std = pvals[(ES > 0 & geZeroMean != 0) | (ES <= 0 & leZeroMean != 0),
                     NES := ES / ifelse(ES > 0, geZeroMean, abs(leZeroMean))],
         pos = pvals[(ES >= 0 & geZeroMean != 0), NES := ES / geZeroMean],
         neg = pvals[(ES <= 0 & leZeroMean != 0), NES := ES / abs(leZeroMean)])
  
  pvals[, pval := as.numeric(NA)]
  pvals[!is.na(NES), pval := pmin((1+nLeEs) / (1 + nLeZero),
                                  (1+nGeEs) / (1 + nGeZero))]
  
  switch(scoreType,
         std = pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)],
         pos = pvals[, nMoreExtreme :=  nGeEs],
         neg = pvals[, nMoreExtreme :=  nLeEs])
  
  pvals[, padj := as.numeric(NA)]
  pvals[!is.na(pval), padj := p.adjust(pval, method = "BH")]
  
  switch(scoreType,
         std = pvals[, nMoreExtreme :=  ifelse(ES > 0, nGeEs, nLeEs)],
         pos = pvals[, nMoreExtreme :=  nGeEs],
         neg = pvals[, nMoreExtreme :=  nLeEs])
  
  pvals[, size := pathwaysSizes[pathway]]
  pvals[, pathway := names(pathwaysFiltered)[pathway]]
  pvals[, leadingEdge := .(leadingEdges)]
  
  return(pvals)
}


## Function randomSelection

randomSelection <- function(universeSize, randomSelectionSize, withCtrl = FALSE, ctrlVariable = NULL, ctrlDistrib = NULL, 
                            printDebug = FALSE){
  res <- c()

  if (printDebug){
    browser()
  }
  
  # if withCtrl = FALSE, just get a list of integer of the required size from the universe
  if (!withCtrl){
    res <- sample.int(universeSize, randomSelectionSize)
  } else if (withCtrl){
    # else, if withCtrl = TRUE, get a list of integer with the same distribution for the control variable as the ctrl distribution
    res <- c()
    # for each bins in the distribution
    for (bin in names(ctrlDistrib)){
      ctrlVariableIdx <- which(ctrlVariable == bin)
      res <- c(res, sample(ctrlVariableIdx, ctrlDistrib[bin]))
    }
  }
  return(res)
}


## Function to prepare the ranked list of genes

prepare_ranked_list <- function(de_results, constraint_data, removeSexChr = TRUE, removePseudotimeGenes = FALSE, 
                                signed = TRUE, typeAnalysis = "withoutCtrl", useEnsemblID = FALSE, printDebug = FALSE){
  if (printDebug){
    print("----> prepare data")
  }
  
  if (removeSexChr){
    # keep only autosomal genes
    de_results <- de_results[which(!de_results$chr %in% c("chrX", "chrY", "chrM")),]
  }
  if (removePseudotimeGenes){
    # load pseudotime data
    load("/Users/benoicla/Desktop/Workspace/HDBR/sexDE_prenatal_brain/data/2-pseudotime/forebrain_pseudotime_analysis.RData")
    # remove genes used in pseudotime inference
    de_results <- de_results[which(!de_results$gene_name %in% rownames(selectedCounts)),]
  }
  # remove the version number from the ensembl id
  de_results$gene_id <- sapply(strsplit(as.character(de_results$gene_id), "\\."), function(x){x[1]})
  
  # add expression bins
  de_results$exprbin <- cut(de_results$AveExpr, 
                            quantile(de_results$AveExpr, seq(0,1,0.1)), 
                            labels=c(1:10), include.lowest=T)
  
  # add constraint bins with LOEUF
  de_results <- merge(de_results, constraint_data, by.x = "gene_id", by.y = "gene_id", all.x = TRUE)
  de_results$constbin <- cut(de_results$LOEUF, 
                             quantile(de_results$LOEUF, seq(0,1,0.1), na.rm=T), 
                             labels=c(1:10), include.lowest=T)
  
  # score used = t-stat
  if (useEnsemblID){
    ranked_de_results <- deframe(de_results[!is.na(de_results$t),c("ensemblID", "t")])
  } else{
    ranked_de_results <- deframe(de_results[!is.na(de_results$t),c("gene_name", "t")])
  }
  # remove duplicated genes
  deduplicatesLines <- !duplicated(names(ranked_de_results))
  ranked_de_results <- ranked_de_results[deduplicatesLines]
  # sort the genes by increasing t-statistic value (if signed) or absolute t-statistic value (if unsigned)
  if (signed){
    idx_sort <- order(ranked_de_results, decreasing = T)
    ranked_de_results <- ranked_de_results[idx_sort]
  } else{
    idx_sort <- order(abs(ranked_de_results), decreasing = T)
    ranked_de_results <- abs(ranked_de_results[idx_sort])
  }
  # sort data frame of results
  de_results <- de_results[idx_sort,]
  
  if (typeAnalysis == "withoutCtrl"){
    ctrlVar <- NULL
  } else if (typeAnalysis == "ctrlExpr"){
    if (useEnsemblID){
      ctrlVar <- de_results$exprbin[which(de_results$ensemblID %in% names(ranked_de_results))]
      names(ctrlVar) <- de_results$ensemblID
    } else{
      ctrlVar <- de_results$exprbin[which(de_results$gene_name %in% names(ranked_de_results))]
      names(ctrlVar) <- de_results$gene_name
    }
  } else if (typeAnalysis == "ctrlConst"){
    if (useEnsemblID){
      ctrlVar <- de_results$constbin[which(de_results$ensemblID %in% names(ranked_de_results))]
      names(ctrlVar) <- de_results$ensemblID
    } else{
      ctrlVar <- de_results$constbin[which(de_results$gene_name %in% names(ranked_de_results))]
      names(ctrlVar) <- de_results$gene_name
    }
  }
  
  return(list("rnk" = ranked_de_results, "de_results" = de_results, "ctrlVar" = ctrlVar))
}

## Function to prepare the ranked list of genes

prepare_ranked_list_bayes <- function(bayes_results, de_results, constraint_data, typeScore, removeSexChr = TRUE, 
                                      signed = TRUE, typeAnalysis = "withoutCtrl", useEnsemblID = FALSE, printDebug = FALSE){
  if (printDebug){
    print("----> prepare data")
    browser()
  }
  
  # remove the version number from the ensembl id
  de_results$ensemblID <- sapply(strsplit(as.character(de_results$gene_id), "\\."), function(x){x[1]})
  # merge the bayes results with the DE results
  bayes_results_merged <- merge(bayes_results, de_results, by = "ensemblID", suffixes = c("", "_DE"))
  
  if (removeSexChr){
    # keep only autosomal genes
    bayes_results_merged <- bayes_results_merged[which(!bayes_results_merged$chr %in% c("chrX", "chrY", "chrM")),]
  }
  
  # add expression bins
  bayes_results_merged$exprbin <- cut(bayes_results_merged$AveExpr, 
                                      quantile(bayes_results_merged$AveExpr, seq(0,1,0.1)), 
                                      labels=c(1:10), include.lowest=T)
  
  # add constraint bins with LOEUF
  bayes_results_merged <- merge(bayes_results_merged, constraint_data, by.x = "ensemblID", by.y = "gene_id", all.x = TRUE)
  bayes_results_merged$constbin <- cut(bayes_results_merged$LOEUF, 
                                       quantile(bayes_results_merged$LOEUF, seq(0,1,0.1), na.rm=T), 
                                       labels=c(1:10), include.lowest=T)
  
  # compute the score
  if (typeScore == "shared"){
    if (signed){
      bayes_results_merged$score <- sign(bayes_results_merged$logFC) * (bayes_results_merged$SHARED0.5 + bayes_results_merged$SHARED1 + bayes_results_merged$SHARED2)
    } else{
      bayes_results_merged$score <- (bayes_results_merged$SHARED0.5 + bayes_results_merged$SHARED1 + bayes_results_merged$SHARED2)
    }
  } else if (typeScore == "prenatal specific"){
    if (signed){
      bayes_results_merged$score <- sign(bayes_results_merged$logFC) * bayes_results_merged$PRENATAL
    } else{
      bayes_results_merged$score <- bayes_results_merged$PRENATAL
    }
  } else if (typeScore == "adult specific"){
    if (signed){
      bayes_results_merged$score <- sign(bayes_results_merged$logFC) * bayes_results_merged$ADULT
    } else{
      bayes_results_merged$score <- bayes_results_merged$ADULT
    }
  }
  
  # score used = bayesian probability
  if (useEnsemblID){
    ranked_de_results <- deframe(bayes_results_merged[!is.na(bayes_results_merged$score), c("ensemblID", "score")])
  } else{
    ranked_de_results <- deframe(bayes_results_merged[!is.na(bayes_results_merged$score), c("gene_name", "score")])
  }
  # remove duplicated genes
  deduplicatesLines <- !duplicated(names(ranked_de_results)) # there is no duplicated lines
  ranked_de_results <- ranked_de_results[deduplicatesLines]
  
  # sort the genes by increasing t-statistic value (if signed) or absolute t-statistic value (if unsigned)
  if (signed){
    idx_sort <- order(ranked_de_results, decreasing = T)
    ranked_de_results <- ranked_de_results[idx_sort]
  } else{
    idx_sort <- order(abs(ranked_de_results), decreasing = T)
    ranked_de_results <- abs(ranked_de_results[idx_sort])
  }
  # sort data frame of results
  bayes_results_merged <- bayes_results_merged[idx_sort,]
  
  if (typeAnalysis == "withoutCtrl"){
    ctrlVar <- NULL
  } else if (typeAnalysis == "ctrlExpr"){
    if (useEnsemblID){
      ctrlVar <- bayes_results_merged$exprbin[which(bayes_results_merged$ensemblID %in% names(ranked_de_results))]
      names(ctrlVar) <- bayes_results_merged$ensemblID
    } else{
      ctrlVar <- bayes_results_merged$exprbin[which(bayes_results_merged$gene_name %in% names(ranked_de_results))]
      names(ctrlVar) <- bayes_results_merged$gene_name
    }
  } else if (typeAnalysis == "ctrlConst"){
    if (useEnsemblID){
      ctrlVar <- bayes_results_merged$constbin[which(bayes_results_merged$ensemblID %in% names(ranked_de_results))]
      names(ctrlVar) <- bayes_results_merged$ensemblID
    } else{
      ctrlVar <- bayes_results_merged$constbin[which(bayes_results_merged$gene_name %in% names(ranked_de_results))]
      names(ctrlVar) <- bayes_results_merged$gene_name
    }
  }
  
  return(list("rnk" = ranked_de_results, "de_results" = bayes_results_merged, "ctrlVar" = ctrlVar))
}

## Function to launch the full enrichment analysis on a one dataset

enrichment1dataset <- function(de_results, typeData, typeAnalysis, genesets, constraint_data, bayes_results,
                               genesetName = "", dataName = "dataset", signed = TRUE, useEnsemblID = FALSE, 
                               outputResults, printDebug = FALSE){
  # init seed
  set.seed(100)
  
  if (printDebug){
    browser()
  }
  
  # prepare data
  if (typeData == "sexDE"){
    preparedData <- prepare_ranked_list(de_results, removeSexChr = TRUE, constraint_data = constraint_data,
                                        removePseudotimeGenes = FALSE, signed = signed, typeAnalysis = typeAnalysis,
                                        useEnsemblID = useEnsemblID, printDebug = printDebug)
    ranked_de_results <- preparedData$rnk
    de_results <- preparedData$de_results
    ctrlVar <- preparedData$ctrlVar
    # colors for lollipop plot
    color <- c("0" = "#919191", "1" = "#074c9b", "2" = "#B2182B")
  } else if (typeData == "pseudotimeDE"){
    preparedData <- prepare_ranked_list(de_results, removeSexChr = TRUE, constraint_data = constraint_data,
                                        removePseudotimeGenes = FALSE, signed = signed, typeAnalysis = typeAnalysis,
                                        useEnsemblID = useEnsemblID, printDebug = printDebug)
    ranked_de_results <- preparedData$rnk
    de_results <- preparedData$de_results
    ctrlVar <- preparedData$ctrlVar
    # colors for lollipop plot
    color <- c("0" = "grey", "1" = "chocolate", "2" = 'aquamarine4')
  } else if (typeData %in% c("shared", "prenatal specific", "adult specific")){
    preparedData <- prepare_ranked_list_bayes(bayes_results, de_results, removeSexChr = TRUE, typeScore = typeData,
                                              constraint_data = constraint_data, signed = signed, typeAnalysis = typeAnalysis,
                                              useEnsemblID = useEnsemblID, printDebug = printDebug)
    ranked_de_results <- preparedData$rnk
    de_results <- preparedData$de_results
    ctrlVar <- preparedData$ctrlVar
    # colors for lollipop plot
    color <- c("0" = "#919191", "1" = "#074c9b", "2" = "#B2182B")
  }
  
  # analysis
  scoreType <- ifelse(signed, "std", "pos")
  if (typeAnalysis == "withoutCtrl"){
    # gseaRes <- gsea_withCtrl(genesets, ranked_de_results, nperm = 10000, withCtrl = FALSE,
    #                          scoreType = scoreType)
    gseaRes <- fgsea(genesets, ranked_de_results, nperm = 10000, scoreType = scoreType)
    gseaResTidy <- gseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
  } else if (typeAnalysis == "ctrlExpr"){
    gseaRes <- gsea_withCtrl(genesets, ranked_de_results, nperm = 10000, withCtrl = TRUE,
                             ctrlVariable = ctrlVar, scoreType = scoreType)
    gseaResTidy <- gseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
  } else if (typeAnalysis == "ctrlConst"){
    gseaRes <- gsea_withCtrl(genesets, ranked_de_results, nperm = 10000, withCtrl = TRUE,
                             ctrlVariable = ctrlVar, scoreType = scoreType)
    gseaResTidy <- gseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
  }
  
  # # plot results
  # plotName <- paste0(outputImages, "gsea_",  typeAnalysis, "_", genesetName, "DiseaseGeneSet_",
  #                    dataName, "_", typeData, "_", ifelse(signed, "signed", "unsigned"), ".png")
  # lollipop_plot(gseaResTidy, inverse = TRUE, padj_threshold = 0.05, fileSave = plotName, signed = signed, colors = color)
  # 
  # save results in a file
  resFileName <- paste0(outputResults, "gsea_",  typeAnalysis, "_", genesetName, "DiseaseGeneSet_",
                        dataName, "_", typeData, "_", ifelse(signed, "signed", "unsigned"), ".txt")
  write_res_gsea(gseaResTidy, resFileName)
  
  # return results
  return(gseaResTidy)
}


#---------------------------------------------#
# Functions to visualize and save the results #
#---------------------------------------------#

## Function to order the pathways for the lollipop plot

order_pathways <- function(data, inverse = FALSE){
  # convert pathways to factor and order them by value of NES
  data$pathway <- factor(data$pathway, ordered = TRUE,
                         levels = data$pathway[order(data$NES)])
  # rename the gene lists
  levels(data$pathway) <- list("Constraint genes" = "constraint_genes",
                               "Highly expressed genes\nin prenatal brain" = "prenatal_high_expressed_genes",
                               "Highly expressed genes\nin adult brain" = "adult_high_expressed_genes", 
                               "Alzheimer" = "alzheimer", 
                               "Amyotrophic lateral\nsclerosis" = "ALS",
                               "Multiple sclerosis" = "MS",
                               "Bipolar disorder" = "BP",
                               "Bipolar I disorder" = "BP1",
                               "Bipolar II disorder" = "BP2",
                               "Bipolar disorder\n(including Schizoaffective)" = "BPSchizo",
                               "Bipolar disorder\nwith psychosis" = "BPwithPsycho",
                               "Bipolar disorder\nwithout psychosis" = "BPwithoutPsycho",
                               "Neurodevelopmental disorders\nwith epilepsy" = "epilepsy_ndd",
                               "Dominant DEE symdrome" = "epilepsy_dee",
                               "Non-acquired focal\nepilepsy" = "NAFE",
                               "Genetic generalized\nepilepsy" = "GGE",
                               "Severe developmental and\nepileptic encephalopathies" = "DEE", 
                               "Epilepsy\n(all types)" = "EPI", 
                               "Schizophrenia" = "scz",
                               "Developmental\nDisorder specific" = "DDspe",
                               "Autism specific" = "ASDspe",
                               "Developmental\nDisorder" = "DD",
                               "Neurodevelopmental\nDisorder" = "NDD",
                               "Autism" = "ASD")
  # # order the disease gene lists by NES values
  # disease <- data[which(data$pathway %in% levels(data$pathway)[6:23]),]
  # diseasePathways <- disease$pathway[order(disease$NES, decreasing = inverse)]
  # # same for control gene lists
  # ctrl <- data[which(data$pathway %in% levels(data$pathway)[1:5]),]
  # ctrlPathways <- ctrl$pathway[order(ctrl$NES, decreasing = inverse, na.last = FALSE)]
  # # reorder the whole pathway list
  # data$pathway <- factor(data$pathway,  ordered = TRUE, 
  #                        levels = c(as.character(ctrlPathways), as.character(diseasePathways)))
  
  return(data)
}


## Function to plot the results as a lollipop plot and save it

lollipop_plot <- function(data, fileSave, inverse = FALSE, padj_threshold = 0.05,
                          colors = c("0" = "#919191", "1" = "#074c9b", "2" = "#B2182B"),
                          signed = TRUE){
  dataOrdered <- order_pathways(data, inverse = inverse)
  # add signif sign
  dataOrdered$signif <- ifelse(dataOrdered$padj < padj_threshold, "*", "")
  # plot enrichments
  p <- ggplot(dataOrdered, aes(x = NES, y = pathway))
  p <- p + geom_segment(aes(x = 0, y = pathway, xend = NES, yend = pathway), color = "grey50")
  p <- p + geom_point(aes(color = as.factor(ifelse(rep(signed, length(NES)), 
                                                   ifelse(NES < 0, 1, 2),
                                                   0))),
                      size = 6)
  p <- p + geom_text(aes(x = NES+sign(NES)*0.1, y = pathway, label = signif), 
                     size = 10)
  p <- p + theme_minimal() + guides(color="none", size="none") + theme(text = element_text(size=15))
  p <- p + scale_color_manual(values = colors)
  p <- p + labs(y = "", x = "Normalized enrichment score")
  p
  # save plot
  ggsave(fileSave, plot = p, width = 20, height = 25, units = "cm")
  
  return(p)
}


## Function to write results in a file

write_res_gsea <- function(resGsea, fileName){
  fwrite(resGsea[, 1:7], file = fileName, quote = FALSE,
         row.names = FALSE, sep = "\t")
}


## Function to plot the enrichment curve

plotEnrichment_manual <- function(pathway, stats, p=1, scoreType = c("std", "pos", "neg"), ticksSize=0.2, addGeneName=FALSE) {
  scoreType <- match.arg(scoreType)
  # this function is mostly a copy of the function "plotEnrichment" from the fgsea package
  rnk <- rank(-stats)
  ord <- order(rnk)
  
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ p)
  statsAdj <- statsAdj / max(abs(statsAdj))
  
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  
  gseaRes <- compute_ES(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE, scoreType = scoreType)
  print(gseaRes)
  
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))
  
  diff <- (max(tops) - min(bottoms)) / 8
  
  # Getting rid of NOTEs
  x=y=NULL
  g <- ggplot(toPlot, aes(x=x, y=y)) +
    geom_point(color="green", size=0.1) +
    geom_hline(yintercept=max(tops), colour="red", linetype="dashed") +
    geom_hline(yintercept=min(bottoms), colour="red", linetype="dashed") +
    geom_hline(yintercept=0, colour="black") +
    geom_line(color="green") + theme_bw() +
    geom_segment(data=data.frame(x=pathway),
                 mapping=aes(x=x, y=-diff/2,
                             xend=x, yend=diff/2),
                 size=ticksSize) +
    theme(panel.border=element_blank(),
          panel.grid.minor=element_blank()) +
    labs(x="rank", y="enrichment score")
  if (addGeneName){
    df_label <- data.frame(x = pathway, y = rep(-diff/2-0.01, length(pathway)), label = names(statsAdj[pathway]))
    g <- g + geom_text_repel(data = df_label, aes(x = x, y = y, label = label), 
                             size = 1, angle = 90, direction = x, force = 1, segment.size = 0.1,
                             min.segment.length = 0, ylim = c(-diff/2-0.05, -diff/2), vjust = 1)
    # g <- g + annotate("text", x = pathway, y = -diff/2-0.01, 
    #                   label = names(statsAdj[pathway]), size = 2, angle = 90,
    #                   hjust = 1, vjust = 1)
    g <- g + ylim(min(min(bottoms), -diff/2-0.05), max(tops))
  }
  
  g
}
