#!/bin/bash 

#SBATCH --job-name=fltrstats
#SBATCH -o fltrVCFstats2-%j.out
#SBATCH -p normal
#SBATCH --time=96:00:00

#This is a script to quantify the number of individuals, contigs, and variants in a filtering directory
#Need to grab and sort the vcf files, then quantify each file, and output results to new file

#to run:                 sbatch fltrVCFstats2.sbatch <threads> <prefix of input/output> <gz>
#to run on vcf files:    sbatch fltrVCFstats2.sbatch 40 SpeciesX.A
#to run on vcf.gz files: sbatch fltrVCFstats2.sbatch 40 SpeciesX.A gz
#both arguments are optional

#enable_lmod
module load parallel
module load R/gcc7/3.5.3

THREADS=$1
PREFIX=$2
if [[ -z $THREADS ]]; then
THREADS=20
fi

echo ""; echo `date` "Running fltrVCFstats2 with the following settings..."
echo "          DIRECTORY: $(pwd)"
echo "          THREADS: $THREADS"
echo "          FILE PREFIX: $PREFIX"

# function to get first word in name of every chromosom/contig
function GET_CHROM_PREFIX(){
	local VcfFileName=$1
	if [[ $3 == "gz" ]]; then
		zcat $VcfFileName | \
		cut -f1 | \
		tail -n1 | \
		sed 's/[-_\.]/\t/g' | \
		cut -f1
	else
		head -n100 $VcfFileName | \
		cut -f1 | \
		tail -n1 | \
		sed 's/[-_\.]/\t/g' | \
		cut -f1
	fi
}
export -f GET_CHROM_PREFIX

#make header row for dat file
HEADER=$(echo -e "File	FilterID	FilterOrder	NumInd	NumContigs	NumMonorphicContigs	NumSNPContigs	NumInsContigs	NumDelContigs	NumSNPs	NumIns	NumDel	NumMissingGeno	NumGenoLess10X	NumGeno10-19X	NumGeno20-49X	NumGeno50-99X	NumGeno100-999X	NumGenoGreater999X" )

#aggregate information for dat file
if [[ $3 == "gz" ]]; then
	echo ""; echo `date` "Getting universal chromosome/contig prefix ..."
	echo "          Querying $VCF_FILE ..."
	CHROM_PREFIX=$(GET_CHROM_PREFIX $(ls -tr ${PREFIX}*vcf.gz | head -n1))
	echo "          Prefix: $CHROM_PREFIX"
	echo "          It is assumed that all values in the CHROM column of vcf begin with $CHROM_PREFIX"
	echo "          If this is not the case, then fltrVCF will not work properly."

	echo ""; echo `date` "Making column 1 of 19 ..."
	Col01=$(ls -tr ${PREFIX}*vcf.gz )
	echo ""; echo `date` "Making column 2 of 19 ..."
	Col02=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "echo -e -n {} | tr '.' '\t' | cut -f5-6 ")	#Filter Order
	echo ""; echo `date` "Making column 3 of 19 ..."
	Col04=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zcat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l ")   #num individuals
	echo ""; echo `date` "Making column 4 of 19 ..."
	Col05=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep \"^$CHROM_PREFIX\" {} | cut -f1 | uniq | wc -l ")   #num contigs
	echo ""; echo `date` "Making column 5 of 19 ..."
	Col05b=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep \"^$CHROM_PREFIX.*TYPE=\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs
	echo ""; echo `date` "Making column 6 of 19 ..."
	Col06=$(paste <(echo $Col05 | tr " " "\n") <(echo $Col05b | tr " " "\n") | awk 'BEGIN { OFS = "\t" } { $3 = $1 - $2} 1 ' | cut -f3 ) # num monomorphic contigs
	echo ""; echo `date` "Making column 7 of 19 ..."
	Col07=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep \"^$CHROM_PREFIX.*TYPE=snp\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w SNPs
	echo ""; echo `date` "Making column 8 of 19 ..."
	Col08=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep \"^$CHROM_PREFIX.*TYPE=ins\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Insertions
	echo ""; echo `date` "Making column 9 of 19 ..."
	Col09=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep \"^$CHROM_PREFIX.*TYPE=del\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Deletions
	echo ""; echo `date` "Making column 10 of 19 ..."
	Col10=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -c \"^$CHROM_PREFIX.*TYPE=snp\" {} ")   #num snps
	echo ""; echo `date` "Making column 11 of 19 ..."
	Col11=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -c \"^$CHROM_PREFIX.*TYPE=ins\" {} ")   #num ins
	echo ""; echo `date` "Making column 12 of 19 ..."
	Col12=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -c \"^$CHROM_PREFIX.*TYPE=del\" {} ")   #num del
	echo ""; echo `date` "Making column 13 of 19 ..."
	Col13=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '\./\.:' {} | wc -l ")   #num missing genotypes
	echo ""; echo `date` "Making column 14 of 19 ..."
	Col14=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9]:' {} | wc -l ")   #num genotypes w less than 10x cvg
	echo ""; echo `date` "Making column 15 of 19 ..."
	Col15=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:1[0-9]:' {} | wc -l ")   #num genotypes w 10-19x cvg
	echo ""; echo `date` "Making column 16 of 19 ..."
	Col16=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l ")   #num genotypes w 20-49x cvg
	echo ""; echo `date` "Making column 17 of 19 ..."
	Col17=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l ")   #num genotypes w 50-99x cvg
	echo ""; echo `date` "Making column 18 of 19 ..."
	Col18=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l ")   #num genotypes w >=100x cvg
	echo ""; echo `date` "Making column 19 of 19 ..."
	Col19=$(ls -tr ${PREFIX}*vcf.gz | parallel --no-notice -j $THREADS -k "zgrep -oh '[01]/[01]:[1-9][0-9][0-9][0-9][0-9]*:' {} | wc -l ")   #num genotypes w >=1000x cvg
else
	echo ""; echo `date` "Getting universal chromosome/contig prefix ..."
	echo "          Querying $VCF_FILE ..."
	CHROM_PREFIX=$(GET_CHROM_PREFIX $(ls -tr ${PREFIX}*vcf | head -n1))
	echo "          Prefix: $CHROM_PREFIX"
	echo "          It is assumed that all values in the CHROM column of vcf begin with $CHROM_PREFIX"
	echo "          If this is not the case, then fltrVCF will not work properly."

	echo ""; echo `date` "Making column 1 of 19 ..."
	Col01=$(ls -tr ${PREFIX}*vcf )
	echo ""; echo `date` "Making column 2 of 19 ..."
	Col02=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "echo -e -n {} | tr '.' '\t' | cut -f5-6 ")	#Filter Order
	echo ""; echo `date` "Making column 3 of 19 ..."
	Col04=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "cat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l ")   #num individuals
	echo ""; echo `date` "Making column 4 of 19 ..."
	Col05=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep \"^$CHROM_PREFIX\" {} | cut -f1 | uniq | wc -l ")   #num contigs
	echo ""; echo `date` "Making column 5 of 19 ..."
	Col05b=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep \"^$CHROM_PREFIX.*TYPE=\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs
	echo ""; echo `date` "Making column 6 of 19 ..."
	Col06=$(paste <(echo $Col05 | tr " " "\n") <(echo $Col05b | tr " " "\n") | awk 'BEGIN { OFS = "\t" } { $3 = $1 - $2} 1 ' | cut -f3 ) # num monomorphic contigs
	echo ""; echo `date` "Making column 7 of 19 ..."
	Col07=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep \"^$CHROM_PREFIX.*TYPE=snp\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w SNPs
	echo ""; echo `date` "Making column 8 of 19 ..."
	Col08=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep \"^$CHROM_PREFIX.*TYPE=ins\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Insertions
	echo ""; echo `date` "Making column 9 of 19 ..."
	Col09=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep \"^$CHROM_PREFIX.*TYPE=del\" {} | cut -f1 | uniq | wc -l ")   #num polymorphic contigs w Deletions
	echo ""; echo `date` "Making column 10 of 19 ..."
	Col10=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -c \"^$CHROM_PREFIX.*TYPE=snp\" {} ")   #num snps
	echo ""; echo `date` "Making column 11 of 19 ..."
	Col11=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -c \"^$CHROM_PREFIX.*TYPE=ins\" {} ")   #num ins
	echo ""; echo `date` "Making column 12 of 19 ..."
	Col12=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -c \"^$CHROM_PREFIX.*TYPE=del\" {} ")   #num del
	echo ""; echo `date` "Making column 13 of 19 ..."
	Col13=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '\./\.:' {} | wc -l ")   #num missing genotypes
	echo ""; echo `date` "Making column 14 of 19 ..."
	Col14=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[1-9]:' {} | wc -l ")   #num genotypes w less than 10x cvg
	echo ""; echo `date` "Making column 15 of 19 ..."
	Col15=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:1[0-9]:' {} | wc -l ")   #num genotypes w 10-19x cvg
	echo ""; echo `date` "Making column 16 of 19 ..."
	Col16=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l ")   #num genotypes w 20-49x cvg
	echo ""; echo `date` "Making column 17 of 19 ..."
	Col17=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l ")   #num genotypes w 50-99x cvg
	echo ""; echo `date` "Making column 18 of 19 ..."
	Col18=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l ")   #num genotypes w >=100x cvg
	echo ""; echo `date` "Making column 19 of 19 ..."
	Col19=$(ls -tr ${PREFIX}*vcf | parallel --no-notice -j $THREADS -k "grep -oh '[01]/[01]:[1-9][0-9][0-9][0-9][0-9]*:' {} | wc -l ")   #num genotypes w >=1000x cvg
fi

#output dat file
echo ""; echo `date` "Assembling tidy data table..."
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
echo ""; echo `date` "Plotting results..."
Rscript plot_fltrVCFstats2.R $PREFIX.fltrStats2.dat $PREFIX.fltrStats2.plots.pdf

echo ""; echo `date` "fltrVCFstats2 complete!"
