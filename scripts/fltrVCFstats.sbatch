#!/bin/bash -l

#SBATCH --job-name=fltrstats
#SBATCH -p main
#SBATCH --time=96:00:00
#SBATCH --ntasks=32

#This is a script to quantify the number of individuals, contigs, and variants in a filtering directory
#Need to grab and sort the vcf files, then quantify each file, and output results to new file

#to run:                 sbatch fltrVCFstats.sbatch <threads> <prefix of input/output> <gz>
#to run on vcf files:    sbatch fltrVCFstats.sbatch 40 SpeciesX.A
#to run on vcf.gz files: sbatch fltrVCFstats.sbatch 40 SpeciesX.A gz
#both arguments are optional

enable_lmod
module load parallel

THREADS=$1
PREFIX=$2
if [[ -z $THREADS ]]; then
	THREADS=20
fi

echo -e "File	NumInd	NumContigs	NumSNPs	NumMissingGeno	NumGenoLess10X	NumGeno10-19X	NumGeno20-49X	NumGeno50-99X	NumGeno100-999X" > $PREFIX.fltrStats.dat

if [[ $3 == "gz" ]]; then
	ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "echo -e -n {}'\t' && \  
		zcat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' &&\   #num individuals
		echo -e -n '\t' && \
		zgrep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' && \   #num contigs
		echo -e -n '\t' && \
		zgrep -c '^dDocent' {} | tr -d '\n' && \   #num snps
		echo -e -n '\t' && \
		zgrep -oh '\./\.:' {} | wc -l | tr -d '\n' && \   #num missing genotypes
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[1-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w less than 10x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:1[0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 10-19x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 20-49x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 50-99x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l && \   #num genotypes w 100-999x cvg
		" >> $PREFIX.fltrStats.dat
else
	ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "echo -e -n {}'\t' && \  
		cat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' &&\   #num individuals
		echo -e -n '\t' && \
		grep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' && \   #num contigs
		echo -e -n '\t' && \
		grep -c '^dDocent' {} | tr -d '\n' && \   #num snps
		echo -e -n '\t' && \
		grep -oh '\./\.:' {} | wc -l | tr -d '\n' && \   #num missing genotypes
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[1-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w less than 10x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:1[0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 10-19x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 20-49x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 50-99x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l && \   #num genotypes w 100-999x cvg
		" >> $PREFIX.fltrStats.dat
fi
