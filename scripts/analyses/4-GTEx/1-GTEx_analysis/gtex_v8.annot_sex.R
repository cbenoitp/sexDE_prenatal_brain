args = commandArgs(trailingOnly=T)

tissue=args[1]
inputdir=args[2]
outdir=args[3]

# load packages
library(data.table)
library(ggplot2)


##### Load data #####

# load annotation
ann <- read.table(sprintf('%s/annotations/gencode.v26.GRCh38.genes.short.txt',inputdir),header=F)
colnames(ann) <- c('chr','gene_id','gene_type','gene_name')
print('annotation loaded')

# load counts
gtex <- read.table(sprintf('%s/counts/%s_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.txt.gz',inputdir,tissue),header=T)
ind_id <- sapply(strsplit(colnames(gtex[-c(1:2)]),"\\."), function(x){o <- unlist(x); paste(o[1],o[2],sep="-")})
samp_id <- gsub("\\.", "-", colnames(gtex[-c(1:2)]))
print('counts loaded')

# reorder annotation to match order of gtex counts
ann <- ann[match(gtex[,1],ann[,2]),]

# load phenotypes
d <- fread(sprintf('%s/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', inputdir), header = TRUE)
print('phenotypes loaded')
d[is.na(d$DTHHRDY),'DTHHRDY'] <- 5

# load samples info
samples <- fread(sprintf('%s/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', inputdir), header = TRUE)
print('samples info loaded')

# load tissues info
tissues <- fread(sprintf('%s/annotations/gtex_tissue_colors.csv', inputdir), header = TRUE, sep = ",")
tissues2 <- tissues[,c("tissue_site_detail", "tissue_site_detail_abbr")]
print('tissues info loaded')

# format phenotypes
d2 <- d[ind_id, on='SUBJID', c('SUBJID','AGE','SEX','DTHHRDY')]
SUBJID <- d2[,SUBJID]
SEX <- d2[,SEX]
DTHHRDY <- as.factor(d2[,DTHHRDY])

# format samples info
samples2 <- samples[samp_id, on='SAMPID', c('SAMPID', 'SMTSD')]
samples2 <- merge(samples2, tissues2, all.x = TRUE, by.x = "SMTSD", by.y = "tissue_site_detail")
TISSUE <- samples2[, tissue_site_detail_abbr]
print('phenotypes formated')

##### REANNOTATE THE SEX OF THE SAMPLES #####

# expression of XIST and of the chrY protein-coding genes
xist <- as.vector(t(gtex[ann$chr == "chrX" & ann$gene_name == "XIST", -c(1:2)]))
chry <- gtex[ann$chr == "chrY" & ann$gene_type == "protein_coding", -c(1:2)]
chry_sum <- as.vector(t(colSums(chry, na.rm = TRUE)))

# function to assign sex using expression of XIST and chrY protein-coding genes as filters.
assignSexFromExpression <- function(xist, chry, thresholdX, thresholdY){
  assigned_sex <- ifelse(chry > thresholdY & xist < thresholdX, 
                         1, ifelse(chry_sum < thresholdY & xist > thresholdX, 2, 0))
  return(assigned_sex)
}

# assign sex using the threshold = 100
samples2$assigned_sex <- assignSexFromExpression(xist, chry_sum, thresholdX = 150, thresholdY = 1000)
print('sex assigned')

# contingency table
table(samples2$assigned_sex, SEX)

# function to plot sex on expression
plotSexOnExpression <- function(xist, chry, sex, thresholdX, thresholdY){
  p <- ggplot()
  p <- p + geom_hline(yintercept = thresholdY)
  p <- p + geom_vline(xintercept = thresholdX)
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

# Plot the chrY expression against xist expression coloring by known sex
p <- plotSexOnExpression(xist, chry_sum, SEX, thresholdX = 150, thresholdY = 1000)
# save plot
namePlot <- sprintf('%s/images/known_sex_GTEx_forebrain.png', outdir)
ggsave(namePlot, plot = p, width = 10, height = 10, units = "cm")

# Plot the chrY expression against xist expression coloring by new assigned sex
p <- plotSexOnExpression(xist, chry_sum, samples2$assigned_sex, thresholdX = 150, thresholdY = 1000)
# save plot
namePlot <- sprintf('%s/images/assigned_sex_GTEx_forebrain.png', outdir)
ggsave(namePlot, plot = p, width = 10, height = 10, units = "cm")
print('plot done')


# format and save the full phenotype table
phenotypes <- data.frame(sampid = samp_id, subjid = SUBJID, sex = SEX, dthhrdy = DTHHRDY, tissue = TISSUE, assigned_sex = samples2$assigned_sex)
write.table(phenotypes, file = sprintf("%s/annotations/GTEx_v8_forebrain_formated_phenotype.txt", outdir),
	col.names = TRUE, row.names = FALSE, quote = FALSE)
print('phenotype file writen')
