#!/bin/bash

#SBATCH --job-name=fltr_3.3_mapJ
#SBATCH --output=output_10.8.J.3.3.out
#SBATCH -p normal
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=rhamner@islander.tamucc.edu
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent
module load perl/5.22.0
module load vcftools/0.1.15
module load rad_haplotyper/1.1.5
module load parallel
module load vcflib/1.0


#files needed in your mapping directory: TotalRawSNPs.vcf, popmap, bam files, reference.fasta
#files needed in your filtering directory: filter_hwe_by_pop.pl, rad_haplotyper.pl
#lines 

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MAPdir="/work/hobi/cbird/opihiMicroGeo2013/mapping_J"
CutoffCode="10.8"
Pops="."   #this is a suffix that you've used to delineate your different popmap files ((*See Below))
DataName="opihiMicroGeo2013.J"
NumProc="20"    #number of processors, cbirdq=40, other hpc nodes = 20, birdlab workstations = either 32 or 40

#*
#the populations in the popmap file must contain at least one, and maybe more individuals for the HWE and Haplotyper scripts to work properly
#Consequently, you may need to adjust the pops in the popmap file.  When you do this, add a suffix on the popmap file as follows:
#original popmap:   popmap.7.10
#modified popmap:   popmap.7.10.Site.LifeStage
#if you want to use the original popmap, then Pops="."

#make a list of contigs for splitting up vcf and parallelization
#grep 'dDocent_Contig' $MAPdir/reference.$CutoffCode.fasta | sed -e 's/>//g' | sort | contigs.txt

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 # echo
 # echo `date` "---------------------------FILTER0001: Isolate SNPs-----------------------------"

# #MNPs are assumed to have been elminated in the FreeBayes Settings, so all we have to do is remove TYPE<>snp
# vcffilter -s -f "TYPE = snp & NUMALT = 1" $MAPdir/TotalRawSNPs.$CutoffCode.vcf > $DataName.$CutoffCode.Filter01.vcf
# snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter01.vcf | wc -l ` 
# contigcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter01.vcf | uniq | wc -l `
# # echo Filter01 $snpcnt > $DataName.$CutoffCode$Pops.output
# echo $snpcnt sites retained and 
# grep ^dDocent_ $DataName.$CutoffCode.Filter01.vcf | cut -f1 | uniq | wc -l
# echo contigs retained


#echo
#echo `date` "---------------------------FILTER0002a: Remove Pools Missing Too Much Data-----------------------------"
# IMPORT VCF into R using vcfR and run assessVCFforFiltering.R to identify pools that are missing too much data and update their names below
# vcfremovesamples $DataName.$CutoffCode.Filter01.vcf noPhix-opihiMicrogeo-Cex-2013-OahuMagicIslandVertical-n7-TAMUCC-008_25-ACTGAT-L003 > $DataName.$CutoffCode.Filter02.vcf

#echo
#echo `date` "---------------------------FILTER0002b: Remove Loci Missing Data for Too Many Pools-----------------------------"
# IMPORT VCF into R using vcfR and run assessVCFforFiltering.R to identify loci that are missing data for too many pools update the LociToRemove.txt

# vcftools --vcf $DataName.$CutoffCode.Filter02.vcf --exclude LociToRemove.txt --out $DataName.$CutoffCode.Filter02b.vcf



# echo
# echo `date` "---------------------------FILTER0003: Min Alternate Allele Count Across All Pools-----------------------------"

# vcffilter -s -f "AC > 5" $DataName.$CutoffCode.Filter02b.vcf > $DataName.$CutoffCode.Filter03.vcf
# snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter03.vcf | wc -l ` 
# contigcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter03.vcf | uniq | wc -l `

 echo
 echo `date` "---------------------------FILTER0004: Filter Low Coverage Allelotypes Per Pool-----------------------------"

# IMPORT VCF into R using vcfR and run assessVCFforFiltering.R to choose min and max DP cutoffs

# -g DP refers to Total read depth at the locus per pool
vcffilter -g "DP > 22" $DataName.$CutoffCode.Filter03.vcf > $DataName.$CutoffCode.Filter04.vcf 
mawk '!/#/' $DataName.$CutoffCode.Filter04.vcf | wc -l


echo echo `date` "---------------------------FILTER0005: Filter High Coverage Allelotypes Per Pool-----------------------------"
# IMPORT VCF into R using vcfR and run assessVCFforFiltering.R to choose min and max DP cutoffs

# -g DP refers to Total read depth at the locus per pool
vcffilter -g "DP < 100" $DataName.$CutoffCode.Filter04.vcf > $DataName.$CutoffCode.Filter05.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter05.vcf | wc -l 


echo
echo `date` "---------------------------FILTER0006: -----------------------------"

#This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)
vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s $DataName.$CutoffCode.Filter05.vcf > $DataName.$CutoffCode.Filter06.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter06.vcf | wc -l


echo
echo `date` "---------------------------FILTER07-----------------------------"

#The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
#loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" $DataName.$CutoffCode.Filter06.vcf > $DataName.$CutoffCode.Filter07.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter07.vcf | wc -l



echo
echo `date` "---------------------------FILTER008-----------------------------"

#first, by removing any locus that has a quality score below 1/4 of the read depth.
vcffilter -f "QUAL / DP > 0.25" $DataName.$CutoffCode.Filter07.vcf > $DataName.$CutoffCode.Filter08.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | wc -l



echo
echo `date` "---------------------------FILTER_LAST: 1 rand snp per contig-----------------------------"

########################################################################################################
#Filter_one_random_snp_per_contig.sh $DataName.$CutoffCode.Filter16 
#Script to take random SNP from every contig in a vcffile

#Calculate number of SNPs
Loci=(`mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | wc -l `)

#Generate list of random numbers
seq 1 500000 | shuf | head -$Loci > nq

#create temporary file that has a random number assigned to each SNP in first column
cat <(mawk '/^#/' $DataName.$CutoffCode.Filter08.vcf) <(paste <(mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | cut -f1-5) nq <(mawk '!/#/' $DataName.$CutoffCode.Filter08.vcf | cut -f7- ) )> $DataName.$CutoffCode.temp

#Get name of VCF file
NAME=$(echo $DataName.$CutoffCode.Filter08.vcf | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 

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



