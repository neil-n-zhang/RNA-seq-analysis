#!/bin/bash

datafolder=~/CA9/raw_data

cd $datafolder

fqs=$(ls)



for fq in $fqs

do
	cd ~/CA9/results
	mkdir $fq

	salmon quant -i /home/neil/bioinfo_tools/Reference/salmon_MusGRCm38CA9_index -l A \
	-1 $datafolder/$fq/${fq}_1.fastq.gz -2 $datafolder/$fq/${fq}_2.fastq.gz \
	-p 8 --validateMappings  -o $fq --seqBias \
 	--useVBOpt \
	--numBootstraps 30 
done
