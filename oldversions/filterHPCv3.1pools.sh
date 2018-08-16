#!/bin/bash

#SBATCH --job-name=fltr_J
#SBATCH --output=output_10.8.J.out
#SBATCH -p normal
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=cbirdtamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent
module load perl/5.22.0
module load vcftools/0.1.15
module load rad_haplotyper/1.1.5
module load parallel
module load vcflib/1.0


#v2.3 decomposes complex variants and removes indels first
#v2.2 updates include generic file names and support for different popmap files

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

 # echo
 # echo `date` "---------------------------FILTER0002: Min Alternate Allele Count Across All Pools-----------------------------"

# vcffilter -s -f "AC > 5" $DataName.$CutoffCode.Filter01.vcf > $DataName.$CutoffCode.Filter02.vcf
# snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter02.vcf | wc -l ` 
# contigcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter02.vcf | uniq | wc -l `
# # echo Filter01 $snpcnt > $DataName.$CutoffCode$Pops.output
# echo $snpcnt sites retained and 
# grep ^dDocent_ $DataName.$CutoffCode.Filter02.vcf | cut -f1 | uniq | wc -l
# echo contigs retained

 echo
 echo `date` "---------------------------FILTER0003: Filter Low Coverage Genotypes-----------------------------"

vcffilter -g "DP > 9" $DataName.$CutoffCode.Filter02.vcf > $DataName.$CutoffCode.Filter03.vcf
snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter03.vcf | wc -l ` 
contigcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter03.vcf | uniq | wc -l `
# echo Filter01 $snpcnt > $DataName.$CutoffCode$Pops.output
echo $snpcnt sites retained and 
grep ^dDocent_ $DataName.$CutoffCode.Filter03.vcf | cut -f1 | uniq | wc -l
echo contigs retained








