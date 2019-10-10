#!/bin/bash

#this script isolates the first 150 bp of a contig in a VCF 

# to run:
# bash limitSites.bash 1 150 ../mkVCF/TotalRawSNPs.15.15.vcf TotalRawSNPs.150.15.15.vcf

startBP=$1
endBP=$2
vcfIN=$3
vcfOUT=$4

cat <(grep -Pv '^dDocent_Contig_[1-9][0-9]*' $vcfIN) \
	<(grep -P '^dDocent_Contig_[1-9][0-9]*' $vcfIN | awk -v bp=$startBP '$2 > bp {print ;}' | awk -v bp=$endBP '$2 < bp {print ;}' ) \
	> $vcfOUT