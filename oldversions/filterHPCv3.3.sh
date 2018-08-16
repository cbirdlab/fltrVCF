#!/bin/bash

#SBATCH --job-name=frOp5.7
#SBATCH --output=output_5.7.N.out
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

MAPdir="/work/hobi/llopezdemesa/opihiSK2014/mapping"
CutoffCode="5.7"
Pops=".Site.LifeStage"   #this is a suffix that you've used to delineate your different popmap files ((*See Below))
DataName="opihiSK2014.N"
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
# echo `date` "---------------------------FILTER0001: DecomposeComplexVars-----------------------------"


# #convert our variant calls to SNPs, This will decompose complex variant calls into phased SNP and INDEL genotypes and keep the INFO flags for loci and genotypes
				# vcfallelicprimitives $MAPdir/TotalRawSNPs.$CutoffCode.vcf --keep-info --keep-geno > $DataName.$CutoffCode.Filter01.vcf

vcffilter -f "AB > 0.25 & AB < 0.75 | AB < 0.01 | AB > 0.99" $MAPdir/TotalRawSNPs.$CutoffCode.vcf | vcfallelicprimitives --keep-info --keep-geno | vcffixup - | vcffilter -f "AC > 2" > $DataName.$CutoffCode.Filter01.vcf
snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter01.vcf | wc -l `
echo Filter01 $snpcnt > $DataName.$CutoffCode$Pops.output
echo $snpcnt sites retained


# echo
# echo `date` "---------------------------FILTER0002: Isolate SNPs-----------------------------"

# #remove indels
vcftools --vcf $DataName.$CutoffCode.Filter01.vcf --remove-indels --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter02
snpcnt=`mawk '!/#/' $DataName.$CutoffCode.Filter02.recode.vcf | wc -l `
echo Filter02 $snpcnt >> $DataName.$CutoffCode$Pops.output


# echo
# echo `date` "---------------------------FILTER0003: min-meanDP-----------------------------"
# #remove loci that have no hope of making it into final dataset
cat $DataName.$CutoffCode.Filter02.recode.vcf | vcffixup - > $DataName.$CutoffCode.Filter02.recode.fixed.vcf
vcftools --vcf $DataName.$CutoffCode.Filter02.recode.fixed.vcf --min-meanDP 5  --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter03

echo
echo `date` "---------------------------FILTER004-----------------------------"

#this line is here to bootstrap filtering from here with a new filter id. if in doubt, comment out
#cp opihiSK2014.L.5.7.Filter03.recode.vcf $DataName.$CutoffCode.Filter03.recode.vcf



#filter raw vcf file with vcf tools

vcftools --vcf $DataName.$CutoffCode.Filter03.recode.vcf --max-missing 0.9 --minQ 30 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter04 


# echo
echo `date` "---------------------------FILTER005-----------------------------"
#set min read depth 
vcftools --vcf $DataName.$CutoffCode.Filter04.recode.vcf --minDP 10 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter05 

#if the line above is disabled then run this line, otherwise disable this line
#mv $DataName.$CutoffCode.Filter04.recode.vcf $DataName.$CutoffCode.Filter05.recode.vcf

echo
echo `date` "---------------------------FILTER006-----------------------------"

#this line is here to bootstrap filtering from here with a new filter id. if in doubt, comment out
#cp opihiSK2014.L.3.10.Filter05.recode.vcf $DataName.$CutoffCode.Filter05.recode.vcf

# get rid of individuals that did not sequence well
vcftools --vcf $DataName.$CutoffCode.Filter05.recode.vcf --missing-indv
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

mawk '$5 > 0.2' $DataName.$CutoffCode.out.imiss | cut -f1 > $DataName.$CutoffCode.lowDP-2.indv

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#list of individuals to remove

echo `date` " Individuals with too much missing data:"
cat $DataName.$CutoffCode.lowDP-2.indv

#remove individuals with low reads
vcftools --vcf $DataName.$CutoffCode.Filter05.recode.vcf --remove $DataName.$CutoffCode.lowDP-2.indv --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter06                                             


echo
echo `date` "---------------------------FILTER007-----------------------------"
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#min mean depth of coverage
#inspect the AF values in the vcf.  This will affect the frequency of rare variants
vcftools --vcf $DataName.$CutoffCode.Filter06.recode.vcf --max-missing 0.9 --maf 0.002 --max-maf 0.998 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter07 --min-meanDP 10


echo
echo `date` "---------------------------FILTER008-----------------------------"

#restrict the data to high percentage of individuals that geneotyped
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#print the popmap file to screen
echo "PopMap File Contents:"
#cat $MAPdir/popmap.$CutoffCode$Pops

#create two lists that have just the individual names for each population
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#need to fix Pop1 and Pop2

#mawk '$2 == "Pop1"' $MAPdir/popmap.$CutoffCode > 1.$CutoffCode.keep && mawk '$2 == "Pop2"' $MAPdir/popmap.$CutoffCode > 2.keep
# mawk '$2 == "noPhiX-opihiSK2014-EastMaui1-1000N-A"' $MAPdir/popmap.$CutoffCode$Pops > 1.$CutoffCode.keep
# mawk '$2 == "noPhiX-opihiSK2014-EastMaui1-1000N-J"' $MAPdir/popmap.$CutoffCode > 2.$CutoffCode.keep

##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# vcftools --vcf $DataName.$CutoffCode.Filter04.recode.vcf --keep 1.$CutoffCode.keep --missing-site --out $DataName.$CutoffCode.1
# vcftools --vcf $DataName.$CutoffCode.Filter04.recode.vcf --keep 2.$CutoffCode.keep --missing-site --out $DataName.$CutoffCode.2

#old script above,new replacement script follows

#get the popnames from the popmap
popnames=`mawk '{print $2}' $MAPdir/popmap.$CutoffCode$Pops | sort | uniq `

#
missingDataByPop() {
	
	mawk -v pop=$1 '$2 == pop' $2 > $4.$3.$1.keep
	vcftools --vcf $4.$3.Filter07.recode.vcf --keep $4.$3.$1.keep --missing-site --out $4.$3.$1
}
export -f missingDataByPop
parallel --no-notice "missingDataByPop {} $MAPdir/popmap.$CutoffCode$Pops $CutoffCode $DataName" ::: ${popnames[*]}

#Print Loci with missing data by pop
#parallel --no-notice "echo {} && cat $DataName.$CutoffCode.{}.lmiss && echo" ::: ${popnames[*]}

# remove loci with more than X proportion of missing data
	#if you have mixed sequence lengths, this will affect if the longer regions are typed
	#cat $DataName.$CutoffCode.*.lmiss | mawk '!/CHR/' | mawk '$6 > 0.6' | cut -f1,2 >> $DataName.$CutoffCode.badloci
	cat $DataName.$CutoffCode.*.lmiss | mawk '!/CHR/' | mawk '$6 > 0.2' | cut -f1,2 | sort | uniq > $DataName.$CutoffCode.badloci

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


vcftools --vcf $DataName.$CutoffCode.Filter07.recode.vcf --exclude-positions $DataName.$CutoffCode.badloci --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter08


echo
echo `date` "---------------------------FILTER009-----------------------------"

vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" $DataName.$CutoffCode.Filter08.recode.vcf > $DataName.$CutoffCode.Filter09.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter09.vcf | wc -l



echo
echo `date` "---------------------------FILTER010-----------------------------"

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)

vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s $DataName.$CutoffCode.Filter09.vcf > $DataName.$CutoffCode.Filter10.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter10.vcf | wc -l

#turn on this line if you turn off the filter above
#cat $DataName.$CutoffCode.Filter09.vcf > $DataName.$CutoffCode.Filter10.vcf

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



echo
echo `date` "---------------------------FILTER011-----------------------------"

#The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
#loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" $DataName.$CutoffCode.Filter10.vcf > $DataName.$CutoffCode.Filter11.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter11.vcf | wc -l



echo
echo `date` "---------------------------FILTER012-----------------------------"

#another filter that can be applied is whether or not their is a discrepancy in the properly paired status  for reads supporting reference or alternate alleles.
#Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
#the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s $DataName.$CutoffCode.Filter11.vcf > $DataName.$CutoffCode.Filter12.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter12.vcf | wc -l


#There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.



echo
echo `date` "---------------------------FILTER0013-----------------------------"

#this is here to bootstrap filtering from here with a new filter id
#cp opihiSK2014.L.25.10.Filter12.vcf $DataName.$CutoffCode.Filter12.vcf

#first, by removing any locus that has a quality score below 1/4 of the read depth.
vcffilter -f "QUAL / DP > 0.25" $DataName.$CutoffCode.Filter12.vcf > $DataName.$CutoffCode.Filter13.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter13.vcf | wc -l



echo
echo `date` "---------------------------FILTER0014-----------------------------"

#second, is a multistep process
#The next step is to create a list of the depth of each locus
cut -f8 $DataName.$CutoffCode.Filter13.vcf | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $DataName.$CutoffCode.Filter13.DEPTH
#Then make a list of quality scores
mawk '!/#/' $DataName.$CutoffCode.Filter13.vcf | cut -f1,2,6 > $DataName.$CutoffCode.Filter13.vcf.loci.qual
#Then calculate mean depth
MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $DataName.$CutoffCode.Filter13.DEPTH )
echo $MeanDepth
MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
echo $MeanDepth2
#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
paste $DataName.$CutoffCode.Filter13.vcf.loci.qual $DataName.$CutoffCode.Filter13.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $DataName.$CutoffCode.Filter13.lowQDloci
cat $DataName.$CutoffCode.Filter13.lowQDloci
#Now we can remove those sites and recalculate the depth across loci with VCFtools
vcftools --vcf $DataName.$CutoffCode.Filter13.vcf --site-depth --exclude-positions $DataName.$CutoffCode.Filter13.lowQDloci --out $DataName.$CutoffCode.Filter13
#Now let’s take VCFtools output and cut it to only the depth scores
cut -f3 $DataName.$CutoffCode.Filter13.ldepth > $DataName.$CutoffCode.Filter13.site.depth
#Now let’s calculate the average depth by dividing the above file by the number of individuals
NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $DataName.$CutoffCode.Filter13.vcf )
mawk '!/D/' $DataName.$CutoffCode.Filter13.site.depth | mawk -v x=$NumInd '{print $1/x}' > $DataName.$CutoffCode.meandepthpersite

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
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
# #Loci that have high mean depth are indicative of either paralogs or multicopy loci. Either way we want to remove them. Here, 
# #I’ll remove all loci above a mean depth of 100, this number should be changed. Now we can combine both filters to produce another VCF file
vcftools --vcf  $DataName.$CutoffCode.Filter13.vcf --recode-INFO-all --out $DataName.$CutoffCode.Filter14 --min-meanDP 20 --max-meanDP 300 --exclude-positions $DataName.$CutoffCode.Filter13.lowQDloci --recode    

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   


echo
echo `date` "---------------------------FILTER0015: HWE-----------------------------"

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

##make sure perl script is in work folder
# Typically, errors would have a low p-value (h setting) and would be present in many populations.
perl filter_hwe_by_pop_HPC.pl -v $DataName.$CutoffCode.Filter14.recode.vcf -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter15.HWE -h 0.001 -d $DataName -co $CutoffCode
mawk '!/#/' $DataName.$CutoffCode.Filter15.HWE.recode.vcf | wc -l

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  

echo
echo `date` "---------------------------FILTER0016: MakeHaplotypes-----------------------------"

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
#note, rad haplotyper skips complex polymorphisms by default
#these are the Puritz Filtering Tutorial default settings
#perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.Filter15.HWE.recode.vcf -x $NumProc -mp 1 -u 20 -ml 4 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter16.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima  

#These are the settings that allowed a legit contig to pass radhapler with a mix of 100 and 150 bp sequences for the maui opihiSK2014 data.  Some of the settings might produce very bad results in another context or dataset
perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.Filter15.HWE.recode.vcf -x $NumProc -e -mp 7 -u 20 -ml 40 -h 300 -d 50 -z 4 -m .5 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter16.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima  

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Why and how many loci were removed by rad_haplotyper
Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
parallel --gnu --null "grep {} $DataName.$CutoffCode.stats.out | cut -f1 > $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
parallel --gnu --null "echo -n {}' removed, ' && wc -l $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
echo file $DataName.$CutoffCode.Filter16.Haplotyped.vcf has 
mawk '!/#/' $DataName.$CutoffCode.Filter16.Haplotyped.vcf | wc -l
echo SNPs
  
echo
echo `date` "---------------------------FILTER0016b: filter radhap rejects-----------------------------"

#this will remove the failures detailed in the previous 3 lines
#Purtiz does this in his filtering tutorial, but I'm thinking that rad_haplotyper does this now. Not sure.
#this takes the stats output from rad_haplotyper and filters the vcf created by the filter_hwe_by_pop script
grep FILTERED $DataName.$CutoffCode.stats.out | mawk '!/Complex/' | cut -f1 > $DataName.$CutoffCode.loci.to.remove
grep -vwf <(cut -f1 $DataName.$CutoffCode.loci.to.remove) $DataName.$CutoffCode.Filter15.HWE.recode.vcf > $DataName.$CutoffCode.Filter16b.vcf
mawk '!/#/' $DataName.$CutoffCode.Filter16b.vcf | wc -l

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
##CEB Need to add a filter here to remove paralogous individuals or to mask genotypes for paralogous individual x genotyp combinations
#radhaplotyper creates a log of errors that might be useful
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
echo
echo `date` "---------------------------FILTER0017: 2 Allelic States-----------------------------"

#this is not in Jon's tutorial.  From Brian Stockwell and Amanda Ackiss
#this removes loci with more than 2 alleles, seems like a good idea, may have already been done above, not sure
vcftools --vcf $DataName.$CutoffCode.Filter16b.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter17
grep -v "^#" $DataName.$CutoffCode.Filter17.recode.vcf | cut -f 1 | uniq | wc -l

#mawk '!/#/' $DataName.$CutoffCode.Filter16.recode.vcf  | wc -l 


  
echo
echo `date` "---------------------------FILTER0018a: 1 rand snp per contig-----------------------------"


########################################################################################################
#Filter_one_random_snp_per_contig.sh $DataName.$CutoffCode.Filter16 
#Script to take random SNP from every contig in a vcffile

#Calculate number of SNPs
Loci=(`mawk '!/#/' $DataName.$CutoffCode.Filter17.recode.vcf | wc -l `)

#Generate list of random numbers
seq 1 500000 | shuf | head -$Loci > nq

#create temporary file that has a random number assigned to each SNP in first column
cat <(mawk '/^#/' $DataName.$CutoffCode.Filter17.recode.vcf) <(paste <(mawk '!/#/' $DataName.$CutoffCode.Filter17.recode.vcf | cut -f1-5) nq <(mawk '!/#/' $DataName.$CutoffCode.Filter17.recode.vcf | cut -f7- ) )> $DataName.$CutoffCode.temp

#Get name of VCF file
NAME=$(echo $DataName.$CutoffCode.Filter18.recode.vcf | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 

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


  
echo
echo `date` "---------------------------FILTER0018b: best snp per contig-----------------------------"

#########################################################################################################
#Script to take the best SNP from every contig in a vcffile

NAME=$(echo $DataName.$CutoffCode.Filter18.recode.vcf | sed -e 's/\.recode.*//g') 

cat $1 | mawk 'BEGIN{last_loc = 0} { 
		if ($1 ~/#/) print $0;
		else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
		else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' > $NAME.bestSNPperLoc.vcf

mawk '!/#/' $NAME.randSNPperLoc.vcf  | wc -l 

echo "Filtered VCF file is saved under name" $NAME.bestSNPperLoc.vcf



