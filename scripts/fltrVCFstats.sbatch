#!/bin/bash -l

#SBATCH --job-name=fltrstats
#SBATCH -p cbirdq
#SBATCH --time=96:00:00

#This is a script to quantify the number of individuals, contigs, and variants in a filtering directory
#Need to grab and sort the vcf files, then quantify each file, and output results to new file

enable_lmod
module load parallel

echo -e "File NumInd NumContigs NumSNPs NumMissingGeno NumGenoLess10X NumGenoLess20X NumGenoLess30X NumGenoLess40X NumGenoLess50X NumGenoLess100X NumGenoLess1000X" > filterStats_2.dat


ls -tr *19*vcf.gz | parallel -k "echo -n {}' ' && \  
	zcat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' &&\   #num individuals
	echo -n ' ' && \
	zgrep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' && \   #num contigs
	echo -n ' ' && \
	zgrep -c '^dDocent' {} | tr -d '\n' && \   #num snps
	echo -n ' ' && \
	zgrep '\./\.:' {} | cut -f 10- | tr '\t' '\n' |  grep '\./\.:' | wc -l | tr -d '\n' && \   #num missing genotypes
	echo -n ' ' && \
	zgrep '[01]/[01]:[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:[1-9]:' | wc -l | tr -d '\n' && \   #num genotypes w less than 10x cvg
	echo -n ' ' && \
	zgrep '[01]/[01]:1[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:1[0-9]:' | wc -l | tr -d '\n' && \   #num genotypes w less than 20x cvg
	echo -n ' ' && \
	zgrep '[01]/[01]:2[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:2[0-9]:' | wc -l | tr -d '\n' && \   #num genotypes w less than 30x cvg
	echo -n ' ' && \
	zgrep '[01]/[01]:3[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:3[0-9]:' | wc -l | tr -d '\n' && \   #num genotypes w less than 40x cvg
	echo -n ' ' && \
	zgrep '[01]/[01]:4[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:4[0-9]:' | wc -l | tr -d '\n' && \   #num genotypes w less than 50x cvg
	echo -n ' ' && \
	zgrep '[01]/[01]:4[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:[5-9][0-9]:' | wc -l | tr -d '\n' && \   #num genotypes w less than 100x cvg
	echo -n ' ' && \
	zgrep '[01]/[01]:4[1-9]:' {} | cut -f 10- | tr '\t' '\n' |  grep '[01]/[01]:[1-9][0-9][0-9]:' | wc -l && \   #num genotypes w less than 1000x cvg
	" >> filterStats_2.dat
