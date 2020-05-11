#!/bin/bash

regex="plot_([1-9]+)_([0-9]+)"
for file in ./stat_analysis_results/*.png
do
	filename=$(basename -- "$file")
	if [[ $filename =~ $regex ]]
	then
		chr="${BASH_REMATCH[1]}"
		locus="${BASH_REMATCH[2]}"

		counts_file="./results_combined/sequences/counts_SNP_chr${chr}.txt"
		cat $counts_file | grep $locus
	fi
done
