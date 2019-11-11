#!/bin/bash 

#SBATCH --job-name=sngltns
#SBATCH -p normal
#SBATCH --time=96:00:00

#This is a script to gather info from SNPs with singleton alternate alleles for downstream determination of whether filters were adequate
#Need to grab and sort the vcf files, then quantify each file, and output results to new file

#to run:                 sbatch singletonStats.sbatch <input vcf> <threads>
#to run on vcf files:    sbatch singletonStats.sbatch SpeciesX.A.vcf 40

module load parallel/gcc7/20190222
#module load R/gcc7/3.5.3

THREADS=$2
VCF=$1

if [[ -z $THREADS ]]; then
	THREADS=20
fi

GetSnpData(){
	i=$1
	VCF=$2
	IndNames=($(grep -m1 '^#CHROM' $VCF | cut -f10- ))
	SnpRecord=$(grep ';AC=1;' $VCF | head -n $i | tail -n 1)
	CONTIG=$(echo $SnpRecord | cut -d" " -f1)
	POS=$(echo $SnpRecord | cut -d" " -f2)
	IndivIndex=$(($(echo $SnpRecord | cut -d" " -f10- | tr -s " " "\n" | cut -d: -f1 | grep -Fn 0/1 | cut -d: -f1) - 1))
	IndivName=$(echo ${IndNames[$IndivIndex]})
	DATA=$(echo $SnpRecord | cut -d" " -f10- | tr -s " " "\n" | grep -Fn 0/1 | tr ":" "\t" | tr "," "\t")
	echo -e "$CONTIG\t$POS\t$DATA\t$IndivName"
}
export -f GetSnpData
NumSingletons=$(grep -c ';AC=1;' $VCF)
parallel --record-env
seq 1 $NumSingletons | parallel --no-notice --env _ -kj $THREADS GetSnpData {} $VCF > Singletons.csv

#gather info on frequency of AC=x
seq 1 600 | parallel --no-notice -kj $THREADS "echo -n AC={}, && grep -c ';AC={};' $VCF" > AC_SNPs.csv

