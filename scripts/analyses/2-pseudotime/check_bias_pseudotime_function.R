# Function to re-do a full pseudotime and DE analysis with a random variable
full_analysis_random <- function(x, phenotypes, gtf.gene.length_df, seedNumber){
  # init seed
  set.seed(seedNumber)
  
  # 1- Pseudotime analysis
  print(">>>>> Pseudotime analysis <<<<<")
  
  # Create a random binary variable
  randomVar <- data.frame("ind" = unique(phenotypes$ind),
                          "randomVar" = sample(c(1,2), size = length(unique(phenotypes$ind)), replace = TRUE, prob = c(0.5, 0.5)))
  # merge this with phenotype
  phenotypes <- left_join(phenotypes, randomVar, by = "ind")
  print("Number of individuals with randomVar = 1:")
  print(sum(phenotypes$randomVar == 1))
  print("Number of individuals with randomVar = 2:")
  print(sum(phenotypes$randomVar == 2))
  
  # Prepare metadata
  RandomVar <- as.factor(phenotypes$randomVar)
  DevStage <- phenotypes$newDevStage
  Individual <- as.factor(phenotypes$ind)
  ForebrainRegion <- phenotypes$newBrainRegion
  Batch <- phenotypes$batch
  # metadata
  metadata <- data.frame(RandomVar = RandomVar, DevStage = DevStage, Individual = Individual, ForebrainRegion = ForebrainRegion, Batch = Batch)
  rownames(metadata) <- phenotypes$experiment_name
  
  # Convert count to TPM
  # Merge gene size with gene metadata
  gene_metadata <- merge(x$genes, gtf.gene.length_df)
  x$genes <- gene_metadata
  # Transform the counts data in log(TPM+1)
  tpm <- counts2tpm(x$counts, x$genes)
  logTpm <- log2(tpm+1)
  rownames(logTpm) <- x$genes$gene_name
  
  # Select the genes to use in the phenopath analysis
  # select higly variable genes
  variance <- rowVars(logTpm)
  ntop <- min(2000, dim(logTpm)[1])
  selectntop <- order(variance, decreasing=TRUE)[seq_len(ntop)]
  selectedCounts <- logTpm[selectntop,]
  # initialize sce object
  sce <- SummarizedExperiment(assays = list(exprs = selectedCounts), 
                              colData = metadata)
  
  # Phenopath analysis
  x_mat <- model.matrix(~ RandomVar)
  fit <- phenopath(sce, x_mat, elbo_tol = 1e-6, thin = 20)
  # Add pseudotime to metadata
  metadata$pseudotime <- trajectory(fit)
  
  # 2- DE analysis
  print(">>>>> DE analysis <<<<<")
  
  # Prepare metadata
  metadata <- metadata[,c(1,6,3)]
  # design for voomDream (the variable to be tested must be a fixed effect)
  form <- ~ RandomVar + pseudotime + RandomVar:pseudotime + (1|Individual)
  # design for sva computation
  design <- model.matrix(~metadata$RandomVar + metadata$pseudotime + metadata$RandomVar:metadata$pseudotime)
  
  # Prepare counts
  suppressWarnings( vobjDream <- voomWithDreamWeights(x, form, metadata, plot = FALSE) ) # Estimate weights using linear mixed model of dream
  # Estimate surrogates variables
  n.sv <- num.sv(vobjDream$E, design, method = "be")
  sva.obj <- sva(vobjDream$E, design, method = "irw", n.sv = n.sv)
  # add svas to metadata
  metadata <- cbind(metadata, sva.obj$sv)
  names(metadata) <- c("RandomVar", "pseudotime", "Individual", paste0("X", 1:n.sv))
  # include svas in formula
  form <- reformulate(c("RandomVar", "pseudotime", "RandomVar:pseudotime", "(1|Individual)", paste0("X", 1:n.sv)))
  # Redo voom normalization with new design matrix (containing svas)
  suppressWarnings( y <- voomWithDreamWeights(x, form, metadata, plot = FALSE) )
  
  # DE analysis
  # Fit the dream model on each gene
  suppressWarnings( fitmm <- dream(y, form, metadata) )
  fitmm2 <- topTable(fitmm, coef = 2, number = 30000, sort.by = "none") 
  fitmm2$qvalue <- qvalue(fitmm2$P.Value)$qvalues
  # compute_mean_expr
  fitmm2$Male <- rowMeans(y$E[,which(metadata$RandomVar==1)])
  fitmm2$VarMale <- rowVars(y$E[,which(metadata$RandomVar==1)])
  fitmm2$Female <- rowMeans(y$E[,which(metadata$RandomVar==2)])
  fitmm2$VarFemale <- rowVars(y$E[,which(metadata$RandomVar==2)])
  
  # Extract DE results
  # RandomVar DE
  de_results <- fitmm2
  # Pseudotime DE
  fitmm_pseudotime <- topTable(fitmm, coef = "pseudotime", number = 30000, sort.by = "none") 
  fitmm_pseudotime$qvalue <- qvalue(fitmm_pseudotime$P.Value)$qvalues
  de_results_pseudotime <- fitmm_pseudotime
  
  # 3- Comparison with main pseudotime-DE results
  # Filter DE results
  genes2remove <- rownames(selectedCounts)
  de_results_pseudotime_filtered <- de_results_pseudotime[which(!de_results_pseudotime$gene_name %in% genes2remove),]
  # Set qvalue threshold
  threshold_qvalue <- 0.01
  
  # Stats
  # Number of DE genes
  nbr_DE_signif <- sum(de_results$qvalue < threshold_qvalue)
  percent_DE_signif <- nbr_DE_signif/length(de_results$gene_name)*100
  print(sprintf("Number of DE genes = %s (representing %.3f %% of the genes)", nbr_DE_signif, percent_DE_signif))
  # down-regulated genes
  nbr_downDE_signif <- length(de_results$gene_name[which(de_results$logFC < 0 & de_results$qvalue < threshold_qvalue)])
  percent_downDE_signif <- nbr_downDE_signif/nbr_DE_signif*100
  print(sprintf("Number of downregulated genes = %s (%.3f %% of DE genes)", nbr_downDE_signif, percent_downDE_signif))
  # upregulated genes
  nbr_upDE_signif <- length(de_results$gene_name[which(de_results$logFC > 0 & de_results$qvalue < threshold_qvalue)])
  percent_upDE_signif <- nbr_upDE_signif/nbr_DE_signif*100
  print(sprintf("Number of upregulated genes = %s (%.3f %% of DE genes)", nbr_upDE_signif, percent_upDE_signif))
  # Number of pseudotime-DE genes
  nbr_sexDE_pseudotime_signif <- sum(de_results_pseudotime_filtered$qvalue < threshold_qvalue)
  percent_sexDE_pseudotime_signif <- nbr_sexDE_pseudotime_signif/length(de_results_pseudotime_filtered$gene_name)*100
  print(sprintf("Number of pseudotime-DE genes = %s (representing %.3f %% of the genes)", nbr_sexDE_pseudotime_signif, percent_sexDE_pseudotime_signif))
  # up-regulated genes
  nbr_up_pseudotime_signif <- length(de_results_pseudotime_filtered$gene_name[which(de_results_pseudotime_filtered$logFC > 0 & de_results_pseudotime_filtered$qvalue < threshold_qvalue)])
  percent_up_pseudotime_signif <- nbr_up_pseudotime_signif/nbr_sexDE_pseudotime_signif*100
  print(sprintf("Number of up-regulated genes = %s (%.3f %% of pseudotime-DE genes)", nbr_up_pseudotime_signif, percent_up_pseudotime_signif))
  # down-regulated genes
  nbr_down_pseudotime_signif <- length(de_results_pseudotime_filtered$gene_name[which(de_results_pseudotime_filtered$logFC < 0 & de_results_pseudotime_filtered$qvalue < threshold_qvalue)])
  percent_down_pseudotime_signif <- nbr_down_pseudotime_signif/nbr_sexDE_pseudotime_signif*100
  print(sprintf("Number of down-regulated genes = %s (%.3f %% of pseudotime-DE genes)", nbr_down_pseudotime_signif, percent_down_pseudotime_signif))
  
  # Compare pseudotime-DE genes from this analysis from the ones observed in the main analysis
  # Load and format pseudotime-DE genes from the main analysis
  de_results_pseudotime_main <- read.table("data/3-DE/VoomDream_topTable_forebrain_pseudotime_sva_ind_interaction_pseudotime.txt",
                                           header = TRUE)
  tmp_env <- new.env()
  load("data/2-pseudotime/forebrain_pseudotime_analysis.RData", envir = tmp_env)
  selectedCounts_main <- tmp_env$selectedCounts
  rm(tmp_env)
  genes2remove <- rownames(selectedCounts_main)
  de_results_pseudotime_main_filtered <- de_results_pseudotime_main[which(!de_results_pseudotime_main$gene_name %in% genes2remove),]
  
  # Overlap of pseudotime-DE genes between the 2 analyses
  # select signif genes
  de_results_pseudotime_filtered_signif <- de_results_pseudotime_filtered[which(de_results_pseudotime_filtered$qvalue < threshold_qvalue), ]
  nbr_pseudoDE_new <- length(de_results_pseudotime_filtered_signif$gene_id)
  print("Number of pseudotime-DE genes in this analysis:")
  print(nbr_pseudoDE_new)
  de_results_pseudotime_main_filtered_signif <- de_results_pseudotime_main_filtered[which(de_results_pseudotime_main_filtered$qvalue < threshold_qvalue), ]
  nbr_pseudoDE_old <- length(de_results_pseudotime_main_filtered_signif$gene_id)
  print("Number of pseudotime-DE genes in main analysis:")
  print(nbr_pseudoDE_old)
  # overlap
  nbr_overlap <- length(intersect(de_results_pseudotime_filtered_signif$gene_id, de_results_pseudotime_main_filtered_signif$gene_id))
  print("Number of common genes:")
  print(nbr_overlap)
  print("Percentage of common genes from main analysis:")
  print(nbr_overlap/nbr_pseudoDE_old*100)
  print("Percentage of common genes from new analysis:")
  print(nbr_overlap/nbr_pseudoDE_new*100)
  nbr_uniq_new <- length(setdiff(de_results_pseudotime_filtered_signif$gene_id, de_results_pseudotime_main_filtered_signif$gene_id))
  print("Number of unique pseudotime-DE genes to the new analysis:")
  print(nbr_uniq_new/nbr_pseudoDE_new*100)
  nbr_uniq_old <- length(setdiff(de_results_pseudotime_main_filtered_signif$gene_id, de_results_pseudotime_filtered_signif$gene_id))
  print("Number of unique pseudotime-DE genes to the main analysis:")
  print(nbr_uniq_old/nbr_pseudoDE_old*100)
  
  # Correlation of logFC
  # merge the 2 pseudotime results
  merged_pseudotime <- merge(de_results_pseudotime[, c(1,6)], de_results_pseudotime_main[, c(1,5)], by = "gene_id", suffixes = c("_new", "_old"))
  # compute correlation
  res_cortest <- cor.test(merged_pseudotime$logFC_new, merged_pseudotime$logFC_old)
  print("Pearson's r =")
  print(res_cortest$estimate)
  print("Pearson's p-value =")
  print(res_cortest$p.value)
  
  return(c("NbrIndRandom1" = sum(phenotypes$randomVar == 1), "NbrIndRandom2" = sum(phenotypes$randomVar == 2), "nbrRandomDE" = nbr_DE_signif,
           "nbrPseudotimeDE" = nbr_pseudoDE_new, "overlap" = nbr_overlap, "percent_overlap" = nbr_overlap/nbr_pseudoDE_old*100,
           "r2" = res_cortest$estimate, "pval" = res_cortest$p.value))
}

# Function to convert count data into log(TPM+1) for pseudotime inference
counts2tpm <- function(counts, genes){
  # RPK = count/geneLength(kb)
  counts <- counts/(genes$size*1e-3)
  # Sum of RPK for each sample
  library_depth <- colSums(counts)
  # TPM = RPK/(sum(RPK)/1e6)
  counts <- sweep(counts, 2, library_depth/1e6,`/`)
  return(counts)
}
