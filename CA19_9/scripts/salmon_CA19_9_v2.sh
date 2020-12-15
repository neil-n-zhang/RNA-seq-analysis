#!/bin/bash

datafolder=$1
outputfolder=$2

cd $datafolder

fqs=$(ls -d */)



for fq in $fqs

do
	fq=$(basename $fq /)
	cd $outputfolder
	mkdir $fq

	salmon quant -i /home/neil/bioinfo_tools/Reference/salmon_MusGRCm38CA9_index -l A \
	-1 $datafolder/$fq/${fq}_1.fastq.gz -2 $datafolder/$fq/${fq}_2.fastq.gz \
	-p 8 --validateMappings  -o $fq --seqBias \
 	--useVBOpt \
	--numBootstraps 30 
done
