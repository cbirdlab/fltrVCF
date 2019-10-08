#!/bin/bash 

#SBATCH --job-name=fltrstats
#SBATCH -p normal
#SBATCH --time=96:00:00

#This is a script to quantify the number of individuals, contigs, and variants in a filtering directory
#Need to grab and sort the vcf files, then quantify each file, and output results to new file

#to run:                 sbatch fltrVCFstats.sbatch <threads> <prefix of input/output> <gz>
#to run on vcf files:    sbatch fltrVCFstats.sbatch 40 SpeciesX.A
#to run on vcf.gz files: sbatch fltrVCFstats.sbatch 40 SpeciesX.A gz
#both arguments are optional

#enable_lmod
module load parallel
module load R/gcc7/3.5.3

THREADS=$1
PREFIX=$2
if [[ -z $THREADS ]]; then
	THREADS=20
fi

#make header row for dat file
HEADER=$(echo -e "File	FilterID	FilterOrder	NumInd	NumContigs	NumSNPs	NumMissingGeno	NumGenoLess10X	NumGeno10-19X	NumGeno20-49X	NumGeno50-99X	NumGeno100-999X" )

#aggregate information for dat file
if [[ $3 == "gz" ]]; then
	Col01=$(ls -tr ${PREFIX}*vcf.gz )
	Col02=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "echo -e -n {} | tr '.' '\t' | cut -f5-6 ")	#Filter Order
	Col04=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zcat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' ")   #num individuals
	Col05=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' ")   #num contigs
	Col06=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -c '^dDocent' {} | tr -d '\n' ")   #num snps
	Col07=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -oh '\./\.:' {} | wc -l | tr -d '\n' ")   #num missing genotypes
	Col08=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w less than 10x cvg
	Col09=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -oh '[01]/[01]:1[0-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w 10-19x cvg
	Col10=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w 20-49x cvg
	Col11=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w 50-99x cvg
	Col12=$(ls -tr ${PREFIX}*vcf.gz | parallel -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l ")   #num genotypes w 100-999x cvg
else
	Col01=$(ls -tr ${PREFIX}*vcf )
	Col02=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "echo -e -n {} | tr '.' '\t' | cut -f5-6 ")	#Filter Order
	Col04=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "cat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' ")   #num individuals
	Col05=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' ")   #num contigs
	Col06=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -c '^dDocent' {} | tr -d '\n' ")   #num snps
	Col07=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -oh '\./\.:' {} | wc -l | tr -d '\n' ")   #num missing genotypes
	Col08=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -oh '[01]/[01]:[1-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w less than 10x cvg
	Col09=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -oh '[01]/[01]:1[0-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w 10-19x cvg
	Col10=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w 20-49x cvg
	Col11=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l | tr -d '\n' ")   #num genotypes w 50-99x cvg
	Col12=$(ls -tr ${PREFIX}*vcf | parallel -j $THREADS -k "grep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l ")   #num genotypes w 100-999x cvg
fi

#output dat file
cat <(echo $HEADER | tr ' ' '\t') \
	<( \
		paste <(echo $Col01 | tr ' ' '\n') \
			<(echo $Col02 | tr ' ' '\n' | paste - -) \
			<(echo $Col04 | tr ' ' '\n') \
			<(echo $Col05 | tr ' ' '\n') \
			<(echo $Col06 | tr ' ' '\n') \
			<(echo $Col07 | tr ' ' '\n') \
			<(echo $Col08 | tr ' ' '\n') \
			<(echo $Col09 | tr ' ' '\n') \
			<(echo $Col10 | tr ' ' '\n') \
			<(echo $Col11 | tr ' ' '\n') \
			<(echo $Col12 | tr ' ' '\n') \
	) > $PREFIX.fltrStats2.dat

#plot dat file
Rscript plot_fltrVCFstats.R $PREFIX.fltrStats2.dat $PREFIX.fltrStats2.pdf
