#!/bin/bash
output_file="mapped_reads.csv"

echo "sample,num_reads" > $output_file

for bam_file in *.bam
do
  sample_name=$(echo $bam_file | cut -d. -f1)
  mapped_reads=$(samtools flagstat $bam_file -@ 8 | grep "with itself and mate mapped" | cut -d " " -f1)

  echo "$sample_name,$mapped_reads" >> $output_file
done
