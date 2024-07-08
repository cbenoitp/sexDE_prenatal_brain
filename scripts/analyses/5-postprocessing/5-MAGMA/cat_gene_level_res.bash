#!/bin/sh

listFiles=`ls /fas/tukiainen/clara/Magma/data/gene_level_analysis/*/*_allLists_MAGMA_res.gsa.out`
finalResFile="/fas/tukiainen/clara/Magma/data/gene_level_analysis/Magma_allLists_res.txt"
for file in $listFiles 
do
	trait=$(basename $(dirname $file))
	echo $trait
	if [ -s "$finalResFile" ]
	then
		#echo "other"
		tail -n 18 $file | sed "s/SET/${trait}/" >> $finalResFile
	else
		#echo "first"
		tail -n 19 $file | sed "s/SET/${trait}/" > $finalResFile
	fi
done
