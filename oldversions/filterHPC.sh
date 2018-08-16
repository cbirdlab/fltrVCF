#!/bin/bash

#SBATCH --job-name=frT119
#SBATCH -p cbirdq
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --mail-user=cbirdtamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load perl/5.22.0
module load vcftools/0.1.15
module load vcflib/1.0
module load rad_haplotyper/1.1.5
module load ddocent/2.2.7
module load parallel


#files needed in your mapping directory: TotalRawSNPs.vcf, popmap, bam files, reference.fasta
#files needed in your filtering directory: filter_hwe_by_pop.pl, rad_haplotyper.py

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MAPdir="/work/GenomicSamples/llopezdemesa/Terrapin2017runs1and2/mapping290"
CutoffCode="11.9"
DataName="Terp290"
NumProc="40"    #number of processors, cbirdq=40, other hpc nodes = 20, birdlab workstations = either 32 or 40

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#filter raw vcf file with vcf tools

vcftools --vcf $MAPdir/TotalRawSNPs.$CutoffCode.vcf --max-missing 0.50 --mac 3 --minQ 30 --recode --recode-INFO-all --out $DataName.$CutoffCode.g5mac3 

#set min depth to 3
vcftools --vcf $DataName.$CutoffCode.g5mac3.recode.vcf --minDP 3 --recode --recode-INFO-all --out $DataName.$CutoffCode.g5mac3dp3

# get rid of individuals that did not sequence well
vcftools --vcf $DataName.$CutoffCode.g5mac3dp3.recode.vcf --missing-indv
rename out. $DataName.$CutoffCode.out. out.imiss
echo; echo "Missing Data Report, Numbers Near 0 are Good"
cat $DataName.$CutoffCode.out.imiss

#graph missing data for individuals
mawk '!/IN/' $DataName.$CutoffCode.out.imiss | cut -f5 > totalmissing
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

mawk '$5 > 0.2' $DataName.$CutoffCode.out.imiss | cut -f1 > $DataName.$CutoffCode.lowDP-2.indv

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#list of individuals to remove

echo `date` " Individuals with too much missing data:"
less $DataName.$CutoffCode.lowDP-2.indv

#remove individuals with low reads
vcftools --vcf $DataName.$CutoffCode.g5mac3dp3.recode.vcf --remove $DataName.$CutoffCode.lowDP-2.indv --recode --recode-INFO-all --out $DataName.$CutoffCode                                             


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#min mean depth of coverage
#restrict the data to high percentage of individuals that geneotyped
vcftools --vcf $DataName.$CutoffCode.recode.vcf --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out $DataName.$CutoffCode.DP3g95maf05 --min-meanDP 10

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#print the popmap file to screen
echo "PopMap File Contents:"
cat $MAPdir/popmap.$CutoffCode

#create two lists that have just the individual names for each population
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#need to fix Pop1 and Pop2
mawk '$2 == "Pop1"' $MAPdir/popmap.$CutoffCode > 1.keep && mawk '$2 == "Pop2"' $MAPdir/popmap.$CutoffCode > 2.keep
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
vcftools --vcf $DataName.$CutoffCode.DP3g95maf05.recode.vcf --keep 1.keep --missing-site --out 1
vcftools --vcf $DataName.$CutoffCode.DP3g95maf05.recode.vcf --keep 2.keep --missing-site --out 2 

rename 1.lmiss $DataName.$CutoffCode.1.lmiss 1.lmiss
rename 2.lmiss $DataName.$CutoffCode.2.lmiss 2.lmiss

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# remove loci with more than X proportion of missing data
cat $DataName.$CutoffCode.1.lmiss $DataName.$CutoffCode.2.lmiss | mawk '!/CHR/' | mawk '$6 > 0.1' | cut -f1,2 >> $DataName.$CutoffCode.badloci

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


vcftools --vcf $DataName.$CutoffCode.DP3g95maf05.recode.vcf --exclude-positions $DataName.$CutoffCode.badloci --recode --recode-INFO-all --out $DataName.$CutoffCode.DP3g95p5maf05

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" $DataName.$CutoffCode.DP3g95p5maf05.recode.vcf > $DataName.$CutoffCode.DP3g95p5maf05.fil1.vcf
mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.recode.vcf | wc -l

mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.fil1.vcf | wc -l


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)

#vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s $DataName.$CutoffCode.DP3g95p5maf05.fil1.vcf > $DataName.$CutoffCode.DP3g95p5maf05.fil2.vcf
#mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.fil2.vcf | wc -l

#turn on this line if you turn off the filter above
cat $DataName.$CutoffCode.DP3g95p5maf05.fil1.vcf > $DataName.$CutoffCode.DP3g95p5maf05.fil2.vcf

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
#loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" $DataName.$CutoffCode.DP3g95p5maf05.fil2.vcf > $DataName.$CutoffCode.DP3g95p5maf05.fil3.vcf
mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.fil3.vcf | wc -l

#another filter that can be applied is whether or not their is a discrepancy in the properly paired status of for reads supporting reference or alternate alleles.
#Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
#the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s $DataName.$CutoffCode.DP3g95p5maf05.fil3.vcf > $DataName.$CutoffCode.DP3g95p5maf05.fil4.vcf
mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.fil4.vcf | wc -l


#There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.

#first, by removing any locus that has a quality score below 1/4 of the read depth.
vcffilter -f "QUAL / DP > 0.25" $DataName.$CutoffCode.DP3g95p5maf05.fil4.vcf > $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf
mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf | wc -l

#second, is a multistep process
#The next step is to create a list of the depth of each locus
cut -f8 $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $DataName.$CutoffCode.DP3g95p5maf05.fil5.DEPTH
#Then make a list of quality scores
mawk '!/#/' $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf | cut -f1,2,6 > $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf.loci.qual
#Then calculate mean depth
MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $DataName.$CutoffCode.DP3g95p5maf05.fil5.DEPTH )
echo $MeanDepth
MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
echo $MeanDepth2
#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
paste $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf.loci.qual $DataName.$CutoffCode.DP3g95p5maf05.fil5.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $DataName.$CutoffCode.DP3g95p5maf05.fil5.lowQDloci
#Now we can remove those sites and recalculate the depth across loci with VCFtools
vcftools --vcf $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf --site-depth --exclude-positions $DataName.$CutoffCode.DP3g95p5maf05.fil5.lowQDloci --out $DataName.$CutoffCode.DP3g95p5maf05.fil5
#Now let’s take VCFtools output and cut it to only the depth scores
cut -f3 $DataName.$CutoffCode.DP3g95p5maf05.fil5.ldepth > $DataName.$CutoffCode.DP3g95p5maf05.fil5.site.depth
#Now let’s calculate the average depth by dividing the above file by the number of individuals
NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf )
mawk '!/D/' $DataName.$CutoffCode.DP3g95p5maf05.fil5.site.depth | mawk -v x=$NumInd '{print $1/x}' > meandepthpersite
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
   
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
#Loci that have high mean depth are indicative of either paralogs or multicopy loci. Either way we want to remove them. Here, 
#I’ll remove all loci above a mean depth of 100, this number should be changed. Now we can combine both filters to produce another VCF file
vcftools --vcf  $DataName.$CutoffCode.DP3g95p5maf05.fil5.vcf --recode-INFO-all --out $DataName.$CutoffCode.DP3g95p5maf05.FIL --max-meanDP 100 --exclude-positions $DataName.$CutoffCode.DP3g95p5maf05.fil5.lowQDloci --recode    

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   
#vcftools --vcf $DataName.$CutoffCode.DP3g95p5maf05.fil5.site.depth --recode-INFO-all --out $DataName.$CutoffCode.vcf   
                                           
###

#HWE Filter
#convert our variant calls to SNPs, This will decompose complex variant calls into phased SNP and INDEL genotypes and keep the INFO flags for loci and genotypes
vcfallelicprimitives $DataName.$CutoffCode.DP3g95p5maf05.FIL.recode.vcf --keep-info --keep-geno > $DataName.$CutoffCode.DP3g95p5maf05.prim.vcf

#remove indels
vcftools --vcf $DataName.$CutoffCode.DP3g95p5maf05.prim.vcf --remove-indels --recode --recode-INFO-all --out $DataName.$CutoffCode.SNP.DP3g95p5maf05



##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##make sure perl script is in work folder
# Typically, errors would have a low p-value (h setting) and would be present in many populations.
perl filter_hwe_by_pop_HPC.pl -v $DataName.$CutoffCode.SNP.DP3g95p5maf05.recode.vcf -p $MAPdir/popmap.$CutoffCode -o $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE -h 0.001 -d $DataName -co $CutoffCode

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.recode.vcf -x $NumProc -mp 1 -u 20 -ml 4 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode -o $DataName.$CutoffCode.SNPS.Haplotyped.vcf -g $DataName.$CutoffCode.haps.genepop -a $DataName.$CutoffCode.haps.ima  

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Why and how many loci were removed by rad_haplotyper
Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
parallel --gnu --null "grep {} $DataName.$CutoffCode.stats.out | cut -f1 > $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
parallel --gnu --null "echo -n {}' removed, ' && wc -l $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"

#this will remove the failures detailed in the previous 3 lines
#Purtiz does this in his filtering tutorial, but I'm thinking that rad_haplotyper does this now. Not sure.
#this takes the stats output from rad_haplotyper and filters the vcf created by the filter_hwe_by_pop script
grep FILTERED $DataName.$CutoffCode.stats.out | mawk '!/Complex/' | cut -f1 > $DataName.$CutoffCode.loci.to.remove
grep -vwf <(cut -f1 $DataName.$CutoffCode.loci.to.remove) $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.recode.vcf > $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.filtered.vcf


#this is not in Jon's tutorial.  From Brian Stockwell and Amanda Ackiss
#this removes loci with more than 2 alleles, seems like a good idea, may have already been done above, not sure
vcftools --vcf $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.filtered.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles
grep -v "^#" $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf | cut -f 1 | uniq | wc -l

#mawk '!/#/' $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf  | wc -l 

########################################################################################################
#Filter_one_random_snp_per_contig.sh $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles 
#Script to take random SNP from every contig in a vcffile

#Calculate number of SNPs
Loci=(`mawk '!/#/' $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf | wc -l `)

#Generate list of random numbers
seq 1 500000 | shuf | head -$Loci > nq

#create temporary file that has a random number assigned to each SNP in first column
cat <(mawk '/^#/' $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf) <(paste <(mawk '!/#/' $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf | cut -f1-5) nq <(mawk '!/#/' $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf | cut -f7- ) )> temp

#Get name of VCF file
NAME=$(echo $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 

#Use awk (mawk) to parse file and select one snp per contig (one with largest random number)
cat temp | mawk 'BEGIN{last_loc = 0} { 
		if ($1 ~/#/) print $0;
		else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' | mawk 'NF > 0' > $NAME.randSNPperLoc.vcf

#Remove temp file
rm temp

#Announce triumphant completion
echo "Filtered VCF file with one random snp per contig is saved as: " $NAME.randSNPperLoc.vcf

#########################################################################################################
#Script to take the best SNP from every contig in a vcffile

NAME=$(echo $DataName.$CutoffCode.SNP.DP3g95p5maf05.HWE.2alleles.recode.vcf | sed -e 's/\.recode.*//g') 

cat $1 | mawk 'BEGIN{last_loc = 0} { 
		if ($1 ~/#/) print $0;
		else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' > $NAME.bestSNPperLoc.vcf


echo "Filtered VCF file is saved under name" $NAME.bestSNPperLoc.vcf



