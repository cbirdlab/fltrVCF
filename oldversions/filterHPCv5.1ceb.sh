#!/bin/bash

#SBATCH --job-name=FilterHPCv5
#SBATCH -p normal
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=tamucc.edu
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent
module load perl/5.22.0
module load vcftools/0.1.15
module load rad_haplotyper/1.1.5
module load parallel
module load vcflib/1.0

#v5.1ceb CEB changed Filter 3 (qual) to be qual per individual, altered filter 6 (Allele Balance) to keep AB = 0, altered filter 7 (min alternate allele count) to also account for all alleles being alternate, removed filter 11 qual/dp>0.25, removed filter 12 qual/dp<2
#v5ceb  CEB reorganized filters and changed settings
#v5 RMH streamlined to match mapping settings in dDocent2.2.4_hpc1+ and dDocent2.2.16_hpc2+
#v4.1 RMH added mechanisms for  "if file exists, skip this step" and a couple on/off switches
#v4 RMH rearranged order 
#v3.2 RMH update to intial decompose complex variants/remove indels that avoids associated genotypes retaining multiple alleles; AB corrected to include 0.99
#v2.3 decomposes complex variants and removes indels first
#v2.2 updates include generic file names and support for different popmap files

# Files needed in your mapping directory:
	# popmap.x.x.xxx, 
	# reference.x.x.fasta
     	# bam files
# Files needed in your filtering directory:
	# filter_hwe_by_pop.pl
	# rad_haplotyperHPC115.pl
        # TotalRawSNPs.x.x.vcf
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MAPdir="/work/vietnam/2018Maymester/cbird/mapping2take6"
CutoffCode="3.6"
Pops=""   #this is a suffix that you've used to delineate your different popmap files ((*See Below))
DataName="Honomanu.F"
NumProc="20"    #number of processors, cbirdq=40, other hpc nodes = 20, birdlab workstations = either 32 or 40

# on/off switches
IndivMissData="on"	# Removing individuals with too much missing data
SitesMissData="on"	# Removing sites with too much missing data
FHWE="on"		# HWE filter
RADHAP="on"		# rad_haplotyper
	
# Popmap file notes
# the populations in the popmap file must contain at least one, and maybe more individuals for the HWE and Haplotyper scripts to work properly
# Consequently, you may need to adjust the pops in the popmap file.  When you do this, add a suffix on the popmap file as follows:
	# original popmap:   popmap.7.10
	# modified popmap:   popmap.7.10.Site.LifeStage
	# if you want to use the original popmap, then Pops=""
# individual names in the popmap file must match file names EXACTLY (everything before the suffix: .F.fq.gz or R2.fq.gz)

#set variables from input
	#number of individuals
	NumInd=$(wc -l ${MAPdir}popmap.${CutoffCode}.${Pops})
	

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


echo
echo `date` "---------------------------FILTER01: Remove Non-Biallelic Sites -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter01.recode.vcf" ]; then
        # remove sites with more or less than 2 alleles and recode info

	# !!!!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$***************
	#  Remember to change back to TotalRawSNPs.$CutoffCode.vcf
	# !!!!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$****************
	
        vcftools --vcf TotalRawSNPs.Chris.$CutoffCode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter01
else
    	echo "Using previously created file"
#    	mawk '!/#/' $DataName.$CutoffCode.Filter01.recode.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER02: Remove Sites with Indels -----------------------------"


if [ ! -e "$DataName.$CutoffCode.Filter02.recode.vcf" ]; then
	# Remove Indels
	vcftools --vcf $DataName.$CutoffCode.Filter01.recode.vcf --remove-indels --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter02
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter02.recode.vcf | wc -l
fi 


echo
echo `date` "---------------------------FILTER03: Remove Sites with Poor Quality -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter03.recode.fixed.vcf" ]; then
	# Remove sites with "Quality value" less than threshold
	vcftools --vcf $DataName.$CutoffCode.Filter02.recode.vcf --minQ 40 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter03
	cat $DataName.$CutoffCode.Filter03.recode.vcf | vcffixup - > $DataName.$CutoffCode.Filter03.recode.fixed.vcf
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter03.recode.fixed.vcf | wc -l
fi 


echo
echo `date` "---------------------------FILTER04: Remove Sites With Low Mean Depth -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter04.recode.vcf" ]; then
	# remove loci with minimum mean depth over all individuals <X
	vcftools --vcf $DataName.$CutoffCode.Filter03.recode.fixed.vcf --min-meanDP 3 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter04
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter04.recode.vcf | wc -l
fi 


echo
echo `date` "---------------------------FILTER05: Remove sites called in <X proportion of individuals -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter05.recode.vcf" ]; then
	# Remove sites called in <X proportion of individuals
	vcftools --vcf $DataName.$CutoffCode.Filter04.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter05
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter05.recode.vcf | wc -l
fi 


echo
echo `date` "---------------------------FILTER06: Remove sites with Average Allele Balance deviating too far from 0.5  -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter06.vcf" ]; then
	vcffilter -s -f "AB > 0.375 & AB < 0.625 | AB = 0" $DataName.$CutoffCode.Filter05.recode.vcf | vcffixup - > $DataName.$CutoffCode.Filter06.vcf
	#snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter06.vcf | wc -l `
	#echo $snpcnt sites retained
	mawk '!/#/' $DataName.$CutoffCode.Filter06.vcf | wc -l
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter06.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER07: Remove sites with Alternate Allele Count <X -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter07.vcf" ]; then
	# remove sites with alternate allele count <2 (across all indivs)
	vcffilter -f "AC > 0 & AN - AC > 0" $DataName.$CutoffCode.Filter06.vcf > $DataName.$CutoffCode.Filter07.vcf
	#snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter07.vcf | wc -l `
	#echo $snpcnt sites retained
	mawk '!/#/' $DataName.$CutoffCode.Filter07.vcf | wc -l
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter07.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER08: Remove sites covered by both F and R reads-----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter08.vcf" ]; then
	# This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)
	vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s $DataName.$CutoffCode.Filter07.vcf > $DataName.$CutoffCode.Filter08.vcf
	mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | wc -l
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER09: Remove sites with low/high ratio of mean mapping quality of alt to ref -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter09.vcf" ]; then
	# The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
	# loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
	vcffilter -f "MQM / MQMR > 0.90 & MQM / MQMR < 1.11" $DataName.$CutoffCode.Filter08.vcf > $DataName.$CutoffCode.Filter09.vcf
	mawk '!/#/' $DataName.$CutoffCode.Filter09.vcf | wc -l
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter09.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER10: Remove sites with one allele only supported by reads not properly paired -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter10.vcf" ]; then
	# another filter that can be applied is whether or not their is a discrepancy in the properly paired status  for reads supporting reference or alternate alleles.
	# Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
	# the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
	vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s $DataName.$CutoffCode.Filter09.vcf > $DataName.$CutoffCode.Filter10.vcf
	mawk '!/#/' $DataName.$CutoffCode.Filter10.vcf | wc -l
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter10.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER11: Remove sites with Quality/DP ratio < 0.25 -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter11.vcf" ]; then
	# There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.
	# first, by removing any locus that has a quality score below 1/4 of the read depth.
	vcffilter -f "QUAL / DP > 0.25" $DataName.$CutoffCode.Filter10.vcf > $DataName.$CutoffCode.Filter11.vcf
	mawk '!/#/' $DataName.$CutoffCode.Filter11.vcf | wc -l
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter11.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER12: Remove sites with Quality > 2xDP -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter12.recode.vcf" ]; then
	# second, is a multistep process
	# create a list of the depth of each locus
	cut -f8 $DataName.$CutoffCode.Filter11.vcf | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $DataName.$CutoffCode.Filter11.DEPTH
	#Then make a list of quality scores
	mawk '!/#/' $DataName.$CutoffCode.Filter11.vcf | cut -f1,2,6 > $DataName.$CutoffCode.Filter11.vcf.loci.qual
	#Then calculate mean depth
	MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $DataName.$CutoffCode.Filter11.DEPTH )
	echo $MeanDepth
	MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
	echo $MeanDepth2
	#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
	paste $DataName.$CutoffCode.Filter11.vcf.loci.qual $DataName.$CutoffCode.Filter11.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $DataName.$CutoffCode.Filter11.lowQDloci
	cat $DataName.$CutoffCode.Filter11.lowQDloci
	#Now we can remove those sites and recalculate the depth across loci with VCFtools
	vcftools --vcf $DataName.$CutoffCode.Filter11.vcf --site-depth --exclude-positions $DataName.$CutoffCode.Filter11.lowQDloci --out $DataName.$CutoffCode.Filter11
	#Now let’s take VCFtools output and cut it to only the depth scores
	cut -f3 $DataName.$CutoffCode.Filter11.ldepth > $DataName.$CutoffCode.Filter11.site.depth
	#Now let’s calculate the average depth by dividing the above file by the number of individuals
	NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $DataName.$CutoffCode.Filter11.vcf )
	mawk '!/D/' $DataName.$CutoffCode.Filter11.site.depth | mawk -v x=$NumInd '{print $1/x}' > $DataName.$CutoffCode.meandepthpersite

	#####################CEB need to figure out how to pass variables into gnuplot
	cp $DataName.$CutoffCode.meandepthpersite meandepthpersite
gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale
	set xrange [10:150] 
	unset label
	set title "Histogram of mean depth per site"
	set ylabel "Number of Occurrences"
	set xlabel "Mean Depth"
	binwidth=1
	bin(x,width)=width*floor(x/width) + binwidth/2.0
	set xtics 5
	plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
	pause -1
EOF
	gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale
	set xrange [160:300] 
	unset label
	set title "Histogram of mean depth per site"
	set ylabel "Number of Occurrences"
	set xlabel "Mean Depth"
	binwidth=1
	bin(x,width)=width*floor(x/width) + binwidth/2.0
	set xtics 5
	plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
	pause -1
EOF
   
	gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale
	set xrange [10:1000] 
	unset label
	set title "Histogram of mean depth per site"
	set ylabel "Number of Occurrences"
	set xlabel "Mean Depth"
	binwidth=1
	bin(x,width)=width*floor(x/width) + binwidth/2.0
	set xtics 5
	plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
	pause -1
EOF
   
	#Loci that have high mean depth are indicative of either paralogs or multicopy loci. Either way we want to remove them. Here, 
	#I’ll remove all loci above a mean depth of 100, this number should be changed. Now we can combine both filters to produce another VCF file
	vcftools --vcf $DataName.$CutoffCode.Filter11.vcf --exclude-positions $DataName.$CutoffCode.Filter11.lowQDloci --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter12 
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter12.recode.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER13: Remove sites with mean DP over all individuals greater than X -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter13.recode.vcf" ]; then
	# Remove sites with mean DP over all individuals greater than X 
	vcftools --vcf $DataName.$CutoffCode.Filter12.recode.vcf --max-meanDP 130 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter13
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter13.recode.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER14: If individual's genotype has DP < X, convert to missing data -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter14.recode.vcf" ]; then
	# Remove genotypes with <X reads
	vcftools --vcf $DataName.$CutoffCode.Filter13.recode.vcf --minDP 10 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter14
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter14.recode.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER15: Remove sites with maf > minor allele frequency > max-maf -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter15.recode.vcf" ]; then
	# Remove sites with minor allele frequency: maf < x < max-maf
		# inspect the AF values in the vcf.  This will affect the frequency of rare variants
	vcftools --vcf $DataName.$CutoffCode.Filter14.recode.vcf --maf 0.005 --max-maf 0.995 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter15
else
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter15.recode.vcf | wc -l	
fi


echo
echo `date` "---------------------------FILTER16: Remove Individuals with Too Much Missing Data-----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter16.recode.vcf" ]; then
# IndivMissData: on/off switch for this filter at beging of script
	if [ $IndivMissData="on" ]; then
		# get rid of individuals that did not sequence well
		vcftools --vcf $DataName.$CutoffCode.Filter15.recode.vcf --missing-indv
		rename out. $DataName.$CutoffCode.out. out.imiss
		echo; echo "Missing Data Report, Numbers Near 0 are Good"
		cat $DataName.$CutoffCode.out.imiss

		#graph missing data for individuals
		mawk '!/IN/' $DataName.$CutoffCode.out.imiss | cut -f5 > $DataName.$CutoffCode.totalmissing
		cp $DataName.$CutoffCode.totalmissing totalmissing
gnuplot << \EOF 
		#filename=system("echo $DataName.$CutoffCode.totalmissing")
		set terminal dumb size 120, 30
		set autoscale 
		unset label
		set title "Histogram of % missing data per individual. Bars to the left are good."
		set ylabel "Number of Occurrences"
		set xlabel "% of missing data"
		#set yr [0:100000]
		binwidth=0.01
		bin(x,width)=width*floor(x/width) + binwidth/2.0
		plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
		pause -1
EOF

		##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#id individuals to remove based upon the proportion of missing data
		#examine the plot and the $DataName.$CutoffCode.out.imiss file to determine if this cutoff is ok
		mawk '$5 > 0.5' $DataName.$CutoffCode.out.imiss | cut -f1 > $DataName.$CutoffCode.lowDP-2.indv
		##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		#list of individuals to remove
		echo `date` " Individuals with too much missing data:"
		cat $DataName.$CutoffCode.lowDP-2.indv

		#remove individuals with low reads
		vcftools --vcf $DataName.$CutoffCode.Filter15.recode.vcf --remove $DataName.$CutoffCode.lowDP-2.indv --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter16                                      
	else
		echo "Filter16 turned off; Filter16 file created is copy of previous filter file"
		cat $DataName.$CutoffCode.Filter15.recode.vcf > $DataName.$CutoffCode.Filter16.recode.vcf
	fi
else 
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter16.recode.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER17: Remove sites with data missing for too many individuals in a population -----------------------------"

if [ ! -e "$DataName.$CutoffCode.Filter17.recode.vcf" ]; then
# restrict the data to loci with a high percentage of individuals that geneotyped
#SitesMissData: on/off switch for this filter at beginning of script
	if [ $SitesMissData="on" ]; then
		# print the popmap file to screen
		echo "PopMap File Contents:"
		cat $MAPdir/popmap.$CutoffCode$Pops

		# get the popnames from the popmap
		popnames=`mawk '{print $2}' $MAPdir/popmap.$CutoffCode$Pops | sort | uniq `

		missingDataByPop() {
			mawk -v pop=$1 '$2 == pop' $2 > $4.$3.$1.keep
			vcftools --vcf $4.$3.Filter16.recode.vcf --keep $4.$3.$1.keep --missing-site --out $4.$3.$1
		}
		export -f missingDataByPop
		parallel --no-notice "missingDataByPop {} $MAPdir/popmap.$CutoffCode$Pops $CutoffCode $DataName" ::: ${popnames[*]}

		# remove loci with more than X proportion of missing data
		# if you have mixed sequence lengths, this will affect if the longer regions are typed
		# UPDATE % missing threshold here	
		cat $DataName.$CutoffCode.*.lmiss | mawk '!/CHR/' | mawk '$6 > 0.5' | cut -f1,2 | sort | uniq > $DataName.$CutoffCode.badloci
		vcftools --vcf $DataName.$CutoffCode.Filter16.recode.vcf --exclude-positions $DataName.$CutoffCode.badloci --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter17
	else
		echo "Filter17 turned off; Filter17 file created is copy of previous filter file"
		cat $DataName.$CutoffCode.Filter16.recode.vcf > $DataName.$CutoffCode.Filter17.recode.vcf
	fi
else 
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter17.recode.vcf | wc -l
fi


echo
echo `date` "---------------------------FILTER18: Remove sites not in HWE (p<0.001) -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter18.HWE.recode.vcf" ]; then
# restrict the data to loci with a high percentage of individuals that geneotyped
# FHWE: on/off switch for this filter at beginning of script 
	if [ $FHWE="on" ]; then
		##make sure filter_hwe_by_pop_HPC.pl script is in working directory
		# Typically, errors would have a low p-value (h setting) and would be present in many populations.
		perl filter_hwe_by_pop_HPC.pl -v $DataName.$CutoffCode.Filter17.recode.vcf -p $MAPdir/popmap.$CutoffCode$Pops -h 0.001 -d $DataName -co $CutoffCode -o $DataName.$CutoffCode.Filter18.HWE
		mawk '!/#/' $DataName.$CutoffCode.Filter18.HWE.recode.vcf | wc -l
	else
		echo "Filter18: HWE turned off; *Filter18.HWE* file created is just a copy of previous filter file"
		cat $DataName.$CutoffCode.Filter17.recode.vcf > $DataName.$CutoffCode.Filter18.HWE.recode.vcf
	fi
else 
	echo "Using previously created file"
#	mawk '!/#/' $DataName.$CutoffCode.Filter18.HWE.recode.vcf | wc -l
fi

echo
echo `date` "---------------------------FILTER19: Make Haplotypes -----------------------------"
if [ ! -e "$DataName.$CutoffCode.Filter19.Haplotyped.vcf" ]; then
# restrict the data to loci with a high percentage of individuals that geneotyped
# RADHAP: on/off switch for this filter at beginning of script
	if [ $RADHAP="on" ]; then
		###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
		#note, rad haplotyper skips complex polymorphisms by default
		#these are the Puritz Filtering Tutorial default settings
		#perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.Filter16.HWE.recode.vcf -x $NumProc -mp 1 -u 20 -ml 4 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter17.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima  
		#RMH settings
		perl rad_haplotyper116HPC.pl -v $DataName.$CutoffCode.Filter18.HWE.recode.vcf -x $NumProc -e -d 50 -mp 10 -u 30 -ml 10 -h 100 -z 0.2 -m 0.6 -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter19.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima

		#Why and how many loci were removed by rad_haplotyper
		Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
		parallel --gnu --null "grep {} $DataName.$CutoffCode.stats.out | cut -f1 > $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
		parallel --gnu --null "echo -n {}' removed, ' && wc -l $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
		echo file $DataName.$CutoffCode.Filter19.Haplotyped.vcf has 
		mawk '!/#/' $DataName.$CutoffCode.Filter19.Haplotyped.vcf | wc -l
		echo SNPs
  	else
		echo "Filter19: rad_haplotyper turned off"
	fi
else 
	echo "Using previously created files"
#	mawk '!/#/' $DataName.$CutoffCode.Filter19.Haplotyped.vcf | wc -l
fi
  
  
echo
echo `date` "---------------------------FILTER20: Select 1 Random SNP per Contig -----------------------------"

#Script to take random SNP from every contig in a vcffile

#Calculate number of SNPs
Loci=(`mawk '!/#/' $DataName.$CutoffCode.Filter19.Haplotyped.vcf | wc -l `)

#Generate list of random numbers
seq 1 500000 | shuf | head -$Loci > nq

#create temporary file that has a random number assigned to each SNP in first column
cat <(mawk '/^#/' $DataName.$CutoffCode.Filter19.Haplotyped.vcf) <(paste <(mawk '!/#/' $DataName.$CutoffCode.Filter19.Haplotyped.vcf | cut -f1-5) nq <(mawk '!/#/' $DataName.$CutoffCode.Filter19.Haplotyped.vcf | cut -f7- ) )> $DataName.$CutoffCode.temp

#Get name of VCF file
NAME=$(echo $DataName.$CutoffCode.Filter19.Haplotyped.vcf | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 

#Use awk (mawk) to parse file and select one snp per contig (one with largest random number)
cat $DataName.$CutoffCode.temp | mawk 'BEGIN{last_loc = 0} { 
		if ($1 ~/#/) print $0;
		else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' | mawk 'NF > 0' > $NAME.randSNPperLoc.vcf

#Remove temp file
rm $DataName.$CutoffCode.temp

mawk '!/#/' $NAME.randSNPperLoc.vcf  | wc -l 

#Announce triumphant completion
echo "Filtered VCF file with one random snp per contig is saved as: " $NAME.randSNPperLoc.vcf

mv nq $DataName.$CutoffCode.nq
mv totalmissing $DataName.$CutoffCode.totalmissing
mv meandepthpersite $DataName.$CutoffCode.meandepthpersite


