#!/bin/bash

#SBATCH --job-name=Filter
#SBATCH --output=output.out
#SBATCH -p cbirdq
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


#v4 RMH rearranged order 
#v3.2 RMH update to intial decompose complex variants/remove indels that avoids associated genotypes retaining multiple alleles; AB corrected to include 0.99
#v2.3 decomposes complex variants and removes indels first
#v2.2 updates include generic file names and support for different popmap files

# Files needed in your mapping directory:
	# TotalRawSNPs.x.x.vcf
	# popmap.x.x.xxx, 
	# reference.x.x.fasta
     	# bam files
# Files needed in your filtering directory:
	# filter_hwe_by_pop.pl
	# rad_haplotyperHPC115.pl

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MAPdir="/YOUR_MAPPING_DIRECTORY_GOES_HERE/"
CutoffCode="3.6"
Pops=".snp2Grp"   #this is a suffix that you've used to delineate your different popmap files ((*See Below))
DataName="DATANAME"
NumProc="40"    #number of processors, cbirdq=40, other hpc nodes = 20, birdlab workstations = either 32 or 40

# Popmap file notes
# the populations in the popmap file must contain at least one, and maybe more individuals for the HWE and Haplotyper scripts to work properly
# Consequently, you may need to adjust the pops in the popmap file.  When you do this, add a suffix on the popmap file as follows:
	# original popmap:   popmap.7.10
	# modified popmap:   popmap.7.10.Site.LifeStage
# if you want to use the original popmap, then Pops=""
# individual names in the popmap file must match file names EXACTLY (everything before the suffix: .F.fq.gz or R2.fq.gz)

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


echo
echo `date` "---------------------------FILTER01: Biallelic only-----------------------------"

# remove sites with more/less than 2 alleles and recode info
vcftools --vcf $MAPdir/TotalRawSNPs.$CutoffCode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter01

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter01.recode.vcf | wc -l `
echo Filter01 $snpcnt


echo
echo `date` "---------------------------FILTER02: Keep only SNPs & MNPs; Convert MNPs to SNPs; remove indels -----------------------------"

# Remove sites that are not SNP or MNP, i.e., 'complex' loci with multiple alleles
vcffilter -s -f "TYPE = snp | TYPE = mnp" $DataName.$CutoffCode.Filter01.recode.vcf > $DataName.$CutoffCode.Filter02a.vcf 

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter02a.vcf | wc -l `
echo $snpcnt sites retained by Filter2a

# Convert MNP calls to SNPs; THIS DOES NOT UPDATE ASSOCIATED INFO COLUMNS (I think this is ok for these)
vcfallelicprimitives --keep-info --keep-geno $DataName.$CutoffCode.Filter02a.vcf > $DataName.$CutoffCode.Filter02b.vcf

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter02b.vcf | wc -l `
echo $snpcnt sites retained by Filter2b

# Remove Indels
vcftools --vcf $DataName.$CutoffCode.Filter02b.vcf --remove-indels --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter02c

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter02c.recode.vcf | wc -l `
echo $snpcnt sites retained by Filter2c


echo
echo `date` "---------------------------FILTER03: Quality -----------------------------"

vcftools --vcf $DataName.$CutoffCode.Filter02c.recode.vcf --minQ 30 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter03

cat $DataName.$CutoffCode.Filter03.recode.vcf | vcffixup - > $DataName.$CutoffCode.Filter03.recode.fixed.vcf

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter03.recode.fixed.vcf | wc -l `
echo $snpcnt sites retained


echo
echo `date` "---------------------------FILTER04: Remove Sites With Low Mean Depth -----------------------------"

# remove loci with minimum mean depth <10
vcftools --vcf $DataName.$CutoffCode.Filter03.recode.fixed.vcf --min-meanDP 10 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter04


echo
echo `date` "---------------------------FILTER05: Remove sites called in <50% individuals -----------------------------"

# Remove sites called in <X% individuals (pools)
vcftools --vcf $DataName.$CutoffCode.Filter04.recode.vcf --max-missing 0.5 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter05


echo
echo `date` "---------------------------FILTER06: Allele Balance, Alternate Allele Count-----------------------------"

# Remove SITES without balance near 0, 0.5, or 1; remove sites with alternate allele count <2 (across all indivs)
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01 | AB > 0.99" $DataName.$CutoffCode.Filter05.recode.vcf | vcffixup - > $DataName.$CutoffCode.Filter06a.vcf

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter06a.vcf | wc -l `
echo $snpcnt sites retained


vcffilter -f "AC > 3" $DataName.$CutoffCode.Filter06a.vcf > $DataName.$CutoffCode.Filter06b.vcf

snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter06b.vcf | wc -l `
echo $snpcnt sites retained


echo
echo `date` "---------------------------FILTER07: Keep only sites covered by F or R reads-----------------------------"

# This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s $DataName.$CutoffCode.Filter06b.vcf > $DataName.$CutoffCode.Filter07.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter07.vcf | wc -l


echo
echo `date` "---------------------------FILTER08: Ratio of mean mapping quality of alt to ref -----------------------------"

# The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
# loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" $DataName.$CutoffCode.Filter07.vcf > $DataName.$CutoffCode.Filter08.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | wc -l


echo
echo `date` "---------------------------FILTER09: Properly paired status -----------------------------"

# another filter that can be applied is whether or not their is a discrepancy in the properly paired status  for reads supporting reference or alternate alleles.
# Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
# the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s $DataName.$CutoffCode.Filter08.vcf > $DataName.$CutoffCode.Filter09.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter09.vcf | wc -l


echo
echo `date` "---------------------------FILTER10: Quality/DP ratio -----------------------------"

# There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.
# first, by removing any locus that has a quality score below 1/4 of the read depth.
vcffilter -f "QUAL / DP > 0.25" $DataName.$CutoffCode.Filter09.vcf > $DataName.$CutoffCode.Filter10.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter10.vcf | wc -l


echo
echo `date` "---------------------------FILTER11: Quality 2xDP -----------------------------"

# second, is a multistep process
# create a list of the depth of each locus
cut -f8 $DataName.$CutoffCode.Filter10.vcf | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $DataName.$CutoffCode.Filter10.DEPTH
#Then make a list of quality scores
mawk '!/#/' $DataName.$CutoffCode.Filter10.vcf | cut -f1,2,6 > $DataName.$CutoffCode.Filter10.vcf.loci.qual
#Then calculate mean depth
MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $DataName.$CutoffCode.Filter10.DEPTH )
echo $MeanDepth
MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
echo $MeanDepth2
#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
paste $DataName.$CutoffCode.Filter10.vcf.loci.qual $DataName.$CutoffCode.Filter10.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $DataName.$CutoffCode.Filter10.lowQDloci
cat $DataName.$CutoffCode.Filter10.lowQDloci
#Now we can remove those sites and recalculate the depth across loci with VCFtools
vcftools --vcf $DataName.$CutoffCode.Filter10.vcf --site-depth --exclude-positions $DataName.$CutoffCode.Filter10.lowQDloci --out $DataName.$CutoffCode.Filter10
#Now let’s take VCFtools output and cut it to only the depth scores
cut -f3 $DataName.$CutoffCode.Filter10.ldepth > $DataName.$CutoffCode.Filter10.site.depth
#Now let’s calculate the average depth by dividing the above file by the number of individuals
NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $DataName.$CutoffCode.Filter10.vcf )
mawk '!/D/' $DataName.$CutoffCode.Filter10.site.depth | mawk -v x=$NumInd '{print $1/x}' > $DataName.$CutoffCode.meandepthpersite

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
vcftools --vcf $DataName.$CutoffCode.Filter10.vcf --max-meanDP 100 --exclude-positions $DataName.$CutoffCode.Filter10.lowQDloci --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter11 


echo
echo `date` "---------------------------FILTER12: minDP (genotypes)-----------------------------"
# Remove genotypes with <X reads
vcftools --vcf $DataName.$CutoffCode.Filter11.recode.vcf --minDP 3 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter12


echo
echo `date` "---------------------------FILTER13: maf, max-maf -----------------------------"

# Remove sites with minor allele frequency: maf < x < max-maf
	# inspect the AF values in the vcf.  This will affect the frequency of rare variants
vcftools --vcf $DataName.$CutoffCode.Filter12.recode.vcf --maf 0.02 --max-maf 0.98 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter13


echo
echo `date` "---------------------------FILTER14: Remove Individuals with Too Much Missing Data-----------------------------"
# get rid of individuals that did not sequence well
vcftools --vcf $DataName.$CutoffCode.Filter13.recode.vcf --missing-indv
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

mawk '$5 > 0.15' $DataName.$CutoffCode.out.imiss | cut -f1 > $DataName.$CutoffCode.lowDP-2.indv

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#list of individuals to remove

echo `date` " Individuals with too much missing data:"
cat $DataName.$CutoffCode.lowDP-2.indv

#remove individuals with low reads
vcftools --vcf $DataName.$CutoffCode.Filter13.recode.vcf --remove $DataName.$CutoffCode.lowDP-2.indv --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter14                                          


echo
echo `date` "---------------------------FILTER15: % missing individuals per population -----------------------------"

# restrict the data to loci with a high percentage of individuals that geneotyped

# print the popmap file to screen
echo "PopMap File Contents:"
cat $MAPdir/popmap.$CutoffCode$Pops

# create lists that have just the individual names for each population
# get the popnames from the popmap
popnames=`mawk '{print $2}' $MAPdir/popmap.$CutoffCode$Pops | sort | uniq `

missingDataByPop() {
	
	mawk -v pop=$1 '$2 == pop' $2 > $4.$3.$1.keep
	vcftools --vcf $4.$3.Filter14.recode.vcf --keep $4.$3.$1.keep --missing-site --out $4.$3.$1
}
export -f missingDataByPop
parallel --no-notice "missingDataByPop {} $MAPdir/popmap.$CutoffCode$Pops $CutoffCode $DataName" ::: ${popnames[*]}

# Print Loci with missing data by pop
# parallel --no-notice "echo {} && cat $DataName.$CutoffCode.{}.lmiss && echo" ::: ${popnames[*]}

# remove loci with more than X proportion of missing data
	# if you have mixed sequence lengths, this will affect if the longer regions are typed
	# update % missing threshold here	
	cat $DataName.$CutoffCode.*.lmiss | mawk '!/CHR/' | mawk '$6 > 0.9' | cut -f1,2 | sort | uniq > $DataName.$CutoffCode.badloci

vcftools --vcf $DataName.$CutoffCode.Filter14.recode.vcf --exclude-positions $DataName.$CutoffCode.badloci --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter15

# If turn off filter16, turn on following line (and vice versa)
### cat $DataName.$CutoffCode.Filter14.recode.vcf > $DataName.$CutoffCode.Filter15.recode.vcf


echo
echo `date` "---------------------------FILTER16: HWE-----------------------------"

##make sure perl script is in work folder
# Typically, errors would have a low p-value (h setting) and would be present in many populations.
perl filter_hwe_by_pop_HPC.pl -v $DataName.$CutoffCode.Filter15.recode.vcf -p $MAPdir/popmap.$CutoffCode$Pops -h 0.001 -d $DataName -co $CutoffCode -o $DataName.$CutoffCode.Filter16.HWE 
mawk '!/#/' $DataName.$CutoffCode.Filter17.HWE.recode.vcf | wc -l


echo
echo `date` "---------------------------FILTER17: MakeHaplotypes-----------------------------"

###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
#note, rad haplotyper skips complex polymorphisms by default
#these are the Puritz Filtering Tutorial default settings
#perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.Filter16.HWE.recode.vcf -x $NumProc -mp 1 -u 20 -ml 4 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter17.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima  

#RMH
perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.Filter16.HWE.recode.vcf -x $NumProc -d 20 -mp 1 -u 20 -ml 4 -h 100 -z 1 -m 0.85 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter17.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima

#Why and how many loci were removed by rad_haplotyper
Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
parallel --gnu --null "grep {} $DataName.$CutoffCode.stats.out | cut -f1 > $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
parallel --gnu --null "echo -n {}' removed, ' && wc -l $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
echo file $DataName.$CutoffCode.Filter17.Haplotyped.vcf has 
mawk '!/#/' $DataName.$CutoffCode.Filter17.Haplotyped.vcf | wc -l
echo SNPs
  
  
echo
echo `date` "---------------------------FILTER18: 1 rand snp per contig-----------------------------"

#Script to take random SNP from every contig in a vcffile

#Calculate number of SNPs
Loci=(`mawk '!/#/' $DataName.$CutoffCode.Filter17.Haplotyped.vcf | wc -l `)

#Generate list of random numbers
seq 1 500000 | shuf | head -$Loci > nq

#create temporary file that has a random number assigned to each SNP in first column
cat <(mawk '/^#/' $DataName.$CutoffCode.Filter17.Haplotyped.vcf) <(paste <(mawk '!/#/' $DataName.$CutoffCode.Filter17.Haplotyped.vcf | cut -f1-5) nq <(mawk '!/#/' $DataName.$CutoffCode.Filter17.Haplotyped.vcf | cut -f7- ) )> $DataName.$CutoffCode.temp

#Get name of VCF file
NAME=$(echo $DataName.$CutoffCode.Filter18.Haplotyped.vcf | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 

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
