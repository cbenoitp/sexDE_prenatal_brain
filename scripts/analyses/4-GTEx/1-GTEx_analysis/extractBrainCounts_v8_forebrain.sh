#!/bin/bash -l

# Sample names 
COLS=`zcat /fas/tukiainen/clara/GTEx/counts/raw_counts/genes/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz | awk 'NR==3 {print;exit}' | tr '\t' '\n'`

# Choose samples from brain tissue (but removing frontal cortex and cerebellar hemisphere samples because they are duplicate region, and keeping only region from the forebrain)
ROWS=$(grep "Brain" /fas/tukiainen/clara/GTEx/annotations/ID_Tissue_v8.txt | grep -E "Amygdala|Anterior cingulate cortex|Caudate|Frontal Cortex|Hippocampus|Hypothalamus|Nucleus accumbens|Putamen" | cut -f1) 

IDX="1,2,"
for K in $ROWS
do
	i=0
	for C in $COLS
	do
		((i++))
		[ $C == $K ] && IDX+="$i,"
	done
done


tissue="brain"
zcat /fas/tukiainen/clara/GTEx/counts/raw_counts/genes/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz | tail -n +3 | cut -f ${IDX%?} > /fas/tukiainen/clara/GTEx/counts/${tissue}_forebrain_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.txt
gzip /fas/tukiainen/clara/GTEx/counts/${tissue}_forebrain_GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.txt

