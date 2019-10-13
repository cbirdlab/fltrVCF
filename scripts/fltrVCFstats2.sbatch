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
HEADER=$(echo -e "File	FilterID	FilterOrder	NumInd	NumContigs	NumMonorphicContigs	NumSNPContigs	NumInsContigs	NumDelContigs	NumSNPs	NumIns	NumDel	NumMissingGeno	NumGenoLess10X	NumGeno10-19X	NumGeno20-49X	NumGeno50-99X	NumGeno100-999X	NumGenoGreater999X" )

#aggregate information for dat file
if [[ $3 == "gz" ]]; then
	Col01=$(ls -tr ${PREFIX}*vcf.gz )
	Col02=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "echo -e -n {} | tr '.' '\t' | cut -f5-6 ")	#Filter Order
	Col04=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zcat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l ")   #num individuals
	Col05=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep '^dDocent' {} | cut -f1 | uniq | wc -l ")   #num contigs
	Col05b=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep '^dDocent.*TYPE=' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs
	Col06=$(paste <(echo $Col05 | tr " " "\n") <(echo $Col05b | tr " " "\n") | awk 'BEGIN { OFS = "\t" } { $3 = $1 - $2} 1 ' | cut -f3 ) # num monomorphic contigs
	Col07=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep '^dDocent.*TYPE=snp' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w SNPs
	Col08=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep '^dDocent.*TYPE=ins' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Insertions
	Col09=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep '^dDocent.*TYPE=del' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Deletions
	Col10=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -c '^dDocent.*TYPE=snp' {} ")   #num snps
	Col11=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -c '^dDocent.*TYPE=ins' {} ")   #num ins
	Col12=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -c '^dDocent.*TYPE=del' {} ")   #num del
	Col13=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '\./\.:' {} | wc -l ")   #num missing genotypes
	Col14=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9]:' {} | wc -l ")   #num genotypes w less than 10x cvg
	Col15=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:1[0-9]:' {} | wc -l ")   #num genotypes w 10-19x cvg
	Col16=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l ")   #num genotypes w 20-49x cvg
	Col17=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l ")   #num genotypes w 50-99x cvg
	Col18=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l ")   #num genotypes w >=100x cvg
	Col19=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9][0-9][0-9][0-9][0-9]*:' {} | wc -l ")   #num genotypes w >=1000x cvg
else
	Col01=$(ls -tr ${PREFIX}*vcf )
	Col02=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "echo -e -n {} | tr '.' '\t' | cut -f5-6 ")	#Filter Order
	Col04=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "cat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l ")   #num individuals
	Col05=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep '^dDocent' {} | cut -f1 | uniq | wc -l ")   #num contigs
	Col05b=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep '^dDocent.*TYPE=' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs
	Col06=$(paste <(echo $Col05 | tr " " "\n") <(echo $Col05b | tr " " "\n") | awk 'BEGIN { OFS = "\t" } { $3 = $1 - $2} 1 ' | cut -f3 ) # num monomorphic contigs
	Col07=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep '^dDocent.*TYPE=snp' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w SNPs
	Col08=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep '^dDocent.*TYPE=ins' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Insertions
	Col09=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep '^dDocent.*TYPE=del' {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Deletions
	Col10=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -c '^dDocent.*TYPE=snp' {} ")   #num snps
	Col11=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -c '^dDocent.*TYPE=ins' {} ")   #num ins
	Col12=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -c '^dDocent.*TYPE=del' {} ")   #num del
	Col13=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '\./\.:' {} | wc -l ")   #num missing genotypes
	Col14=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[1-9]:' {} | wc -l ")   #num genotypes w less than 10x cvg
	Col15=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:1[0-9]:' {} | wc -l ")   #num genotypes w 10-19x cvg
	Col16=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l ")   #num genotypes w 20-49x cvg
	Col17=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l ")   #num genotypes w 50-99x cvg
	Col18=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l ")   #num genotypes w >=100x cvg
	Col19=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[1-9][0-9][0-9][0-9][0-9]*:' {} | wc -l ")   #num genotypes w >=1000x cvg
fi

#output dat file
cat <(echo $HEADER | tr ' ' '\t') \
	<( \
		paste <(echo $Col01 | tr ' ' '\n') \
			<(echo $Col02 | tr ' ' '\n' | paste - - ) \
			<(echo $Col04 | tr ' ' '\n') \
			<(echo $Col05 | tr ' ' '\n') \
			<(echo $Col06 | tr ' ' '\n') \
			<(echo $Col07 | tr ' ' '\n') \
			<(echo $Col08 | tr ' ' '\n') \
			<(echo $Col09 | tr ' ' '\n') \
			<(echo $Col10 | tr ' ' '\n') \
			<(echo $Col11 | tr ' ' '\n') \
			<(echo $Col12 | tr ' ' '\n') \
			<(echo $Col13 | tr ' ' '\n') \
			<(echo $Col14 | tr ' ' '\n') \
			<(echo $Col15 | tr ' ' '\n') \
			<(echo $Col16 | tr ' ' '\n') \
			<(echo $Col17 | tr ' ' '\n') \
			<(echo $Col18 | tr ' ' '\n') \
			<(echo $Col19 | tr ' ' '\n') \
	) > $PREFIX.fltrStats2.dat

#plot dat file
Rscript plot_fltrVCFstats2.R $PREFIX.fltrStats2.dat $PREFIX.fltrStats2.plots.pdf
