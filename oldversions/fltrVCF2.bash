#!/bin/bash

# Files needed:
	# popmap.x.x.xxx, 
	# reference.x.x.fasta
    # bam files
	# vcf file
# Scripts needed in your filtering directory:
	# filter_hwe_by_pop.pl
	# rad_haplotyperHPC115.pl
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo ""

###################################################################################################################
#Define main function
###################################################################################################################

function MAIN(){
	#read in variables
	name=$1[@]
	FILTERS=("${!name}")
	MODE=$2
	CutoffCode=$3
	BAM_PATH=$4
	VCF_FILE=$5
	REF_FILE=$6
	PopMap=$7
	CONFIG_FILE=$8
	HWE_SCRIPT=$9
	RADHAP_SCRIPT=${10}
	DataName=${11}
	NumProc=${12}

	
	
	for i in ${FILTERS[@]}; do
		FILTER $i $MODE $CutoffCode $BAM_PATH $VCF_FILE $REF_FILE $PopMap $CONFIG_FILE $HWE_SCRIPT $RADHAP_SCRIPT $DataName $NumProc
		VCF_FILE=$(ls -t *vcf | head -n 1)
	done
}


###################################################################################################################
#Define filters 
###################################################################################################################

function FILTER(){

	#read in variables
	FILTER_ID=$1
	MODE=$2
	CutoffCode=$3
	BAM_PATH=$4
	VCF_FILE=$5
	REF_FILE=$6
	PopMap=$7
	CONFIG_FILE=$8
	HWE_SCRIPT=$9
	RADHAP_SCRIPT=${10}
	DataName=${11}
	NumProc=${12}

	# on/off switches
	IndivMissData="on"	# Removing individuals with too much missing data
	SitesMissData="on"	# Removing sites with too much missing data
	FHWE="on"		# HWE filter
	RADHAP="on"		# rad_haplotyper
	

	if [ $FILTER_ID == "01" ]; then
		echo; echo `date` "---------------------------FILTER01: Remove sites with Y < alleles < X -----------------------------"
		#get settings from config file
			MIN=($(grep -P '^\t* *01\t* *vcftools\t* *--min-alleles' ${CONFIG_FILE} | sed 's/\t* *01\t* *vcftools\t* *--min-alleles\t* *//g' | sed 's/\t* *#.*//g' )) 
			if [ -z "$MIN" ]; then MIN=2; fi
			MAX=($(grep -P '^\t* *01\t* *vcftools\t* *--max-alleles' ${CONFIG_FILE} | sed 's/\t* *01\t* *vcftools\t* *--max-alleles\t* *//g' | sed 's/\t* *#.*//g' )) 
			if [ -z "$MAX" ]; then MAX=2; fi
		vcftools --vcf $VCF_FILE --min-alleles $MIN --max-alleles $MAX --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter01
	
	elif [ $FILTER_ID == "02" ]; then
		echo; echo `date` "---------------------------FILTER02: Remove Sites with Indels -----------------------------"
		vcftools --vcf $VCF_FILE --remove-indels --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter02
	
	elif [ $FILTER_ID == "03" ]; then
		echo; echo `date` "---------------------------FILTER03: Remove Sites with QUAL < minQ -----------------------------"
		minQ=($(grep -P '^\t* *03\t* *vcftools\t* *--minQ' ${CONFIG_FILE} | sed 's/\t* *03\t* *vcftools\t* *--minQ\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${minQ}" ]; then ${minQ}=30; fi
		vcftools --vcf $VCF_FILE --minQ ${minQ} --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter03
		#cat $DataName.$CutoffCode.${MODE}.Filter03.recode.vcf | vcffixup - > $DataName.$CutoffCode.${MODE}.Filter03.recode.fixed.vcf
		
	elif [ $FILTER_ID == "04" ]; then
		echo; echo `date` "---------------------------FILTER04: Remove Sites With Mean Depth of Coverage < min-meanDP -----------------------------"
		minMeanDP=($(grep -P '^\t* *04\t* *vcftools\t* *--min-meanDP' ${CONFIG_FILE} | sed 's/\t* *04\t* *vcftools\t* *--min-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${minMeanDP}" ]; then ${minMeanDP}=2; fi
		vcftools --vcf $VCF_FILE --min-meanDP ${minMeanDP} --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter04

	elif [ $FILTER_ID == "05" ]; then
		echo; echo `date` "---------------------------FILTER05: Remove sites called in <X proportion of individuals -----------------------------"
		maxMissing=($(grep -P '^\t* *05\t* *vcftools\t* *--max-missing' ${CONFIG_FILE} | sed 's/\t* *05\t* *vcftools\t* *--max-missing\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${maxMissing}" ]; then ${maxMissing}=0.5; fi
		vcftools --vcf $DataName.$CutoffCode.Filter04.recode.vcf --max-missing ${maxMissing} --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter05
	
	elif [ $FILTER_ID == "06" ]; then
		echo; echo `date` "---------------------------FILTER06: Remove sites with Average Allele Balance deviating too far from 0.5 while keeping those with AB=0  -----------------------------"
		minAB=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *min' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${minAB}" ]; then ${minAB}=0.375; fi
		maxAB=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *max' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *max\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${maxAB}" ]; then ${maxAB}=0.625; fi
		vcffilter -s -f "AB > ${minAB} & AB < ${maxAB} | AB = 0" $VCF_FILE | vcffixup - > $DataName.$CutoffCode.${MODE}.Filter06.vcf
		mawk '!/#/' $DataName.$CutoffCode.${MODE}.Filter06.vcf | wc -l
	
	elif [ $FILTER_ID == "07" ]; then
		echo; echo `date` "---------------------------FILTER07: Remove sites with Alternate Allele Count <=X -----------------------------"
		minAC=($(grep -P '^\t* *07\t* *vcffilter\t* *AC\t* *min' ${CONFIG_FILE} | sed 's/\t* *07\t* *vcffilter\t* *AC\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${minAC}" ]; then ${minAC}=1; fi
		vcffilter -f "AC > ${minAC} & AN - AC > ${minAC}" $VCF_FILE > $DataName.$CutoffCode.${MODE}.Filter07.vcf
		mawk '!/#/' $DataName.$CutoffCode.${MODE}.Filter07.vcf | wc -l
		
	elif [ $FILTER_ID == "08" ]; then
		echo; echo `date` "---------------------------FILTER08: Remove sites covered by both F and R reads-----------------------------"
		# This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)
		THRESHOLD=($(grep -P '^\t* *08\t* *vcffilter\t* *SAF\/SAR\t* *min' ${CONFIG_FILE} | sed 's/\t* *08\t* *vcffilter\t* *SAF\/SAR\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=1; fi
		vcffilter -f "SAF / SAR > ${THRESHOLD} & SRF / SRR > ${THRESHOLD} | SAR / SAF > ${THRESHOLD} & SRR / SRF > ${THRESHOLD}" -s $VCF_FILE > $DataName.$CutoffCode.${MODE}.Filter08.vcf
		mawk '!/#/' $DataName.$CutoffCode.${MODE}.Filter08.vcf | wc -l

	elif [ $FILTER_ID == "09" ]; then
		echo; echo `date` "---------------------------FILTER09: Remove sites with low/high ratio of mean mapping quality of alt to ref -----------------------------"
		# The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
		# loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
		THRESHOLD=($(grep -P '^\t* *09\t* *vcffilter\t* *MQM\/MQMR\t* *min' ${CONFIG_FILE} | sed 's/\t* *09\t* *vcffilter\t* *MQM\/MQMR\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.1; fi
		THRESHOLDmin=$(echo 1-${THRESHOLD} | R --vanilla --quiet | sed -n '2s/.* //p')
		THRESHOLDmax=$(echo 1/${THRESHOLDmin} | R --vanilla --quiet | sed -n '2s/.* //p')
		vcffilter -f "MQM / MQMR > 1 - ${THRESHOLDmin} & MQM / MQMR < ${THRESHOLDmax}" $VCF_FILE > $DataName.$CutoffCode.${MODE}.Filter09.vcf
		mawk '!/#/' $DataName.$CutoffCode.${MODE}.Filter09.vcf | wc -l

	elif [ $FILTER_ID == "10" ]; then
		echo; echo `date` "---------------------------FILTER10: Remove sites with one allele only supported by reads not properly paired -----------------------------"
		# another filter that can be applied is whether or not their is a discrepancy in the properly paired status  for reads supporting reference or alternate alleles.
		# Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
		# the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
		vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s $VCF_FILE > $DataName.$CutoffCode.${MODE}.Filter10.vcf
		mawk '!/#/' $DataName.$CutoffCode.${MODE}.Filter10.vcf | wc -l
	
	elif [ $FILTER_ID == "11" ]; then
		echo; echo `date` "---------------------------FILTER11: Remove sites with Quality/DP ratio < 0.25 -----------------------------"
		# There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.
		# first, by removing any locus that has a quality score below 1/4 of the read depth.
		THRESHOLD=($(grep -P '^\t* *11\t* *vcffilter\t* *QUAL\/DP\t* *min' ${CONFIG_FILE} | sed 's/\t* *11\t* *vcffilter\t* *QUAL\/DP\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.25; fi
		vcffilter -f "QUAL / DP > ${THRESHOLD}" $VCF_FILE > $DataName.$CutoffCode.${MODE}.Filter11.vcf
		mawk '!/#/' $DataName.$CutoffCode.${MODE}.Filter11.vcf | wc -l

	elif [ $FILTER_ID == "12" ]; then
		echo; echo `date` "---------------------------FILTER12: Remove sites with Quality > 2xDP -----------------------------"
		# second, is a multistep process
#		THRESHOLD=($(grep -P '^\t* *12\t* *vcffilter\t* *QUAL\/DP max' ${CONFIG_FILE} | sed 's/\t* *12\t* *vcffilter\t* *QUAL\/DP max\t* *//g' | sed 's/\t* *#.*//g' )) 
#		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=2; fi
		# create a list of the depth of each locus
		cut -f8 $VCF_FILE | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $DataName.$CutoffCode.${MODE}.Filter12.DEPTH
		#Then make a list of quality scores
		mawk '!/#/' $VCF_FILE | cut -f1,2,6 > $DataName.$CutoffCode.${MODE}.Filter12.vcf.loci.qual
		#Then calculate mean depth
		MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $DataName.$CutoffCode.${MODE}.Filter12.DEPTH )
		#echo $MeanDepth
		MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
		#echo $MeanDepth2
		#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
		paste $DataName.$CutoffCode.${MODE}.Filter12.vcf.loci.qual $DataName.$CutoffCode.${MODE}.Filter12.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $DataName.$CutoffCode.${MODE}.Filter12.lowQDloci
		#cat $DataName.$CutoffCode.${MODE}.Filter11.lowQDloci
		vcftools --vcf $VCF_FILE --exclude-positions $DataName.$CutoffCode.${MODE}.Filter12.lowQDloci --recode --recode-INFO-all --out $DataName.$CutoffCode.Filter12 

	elif [ $FILTER_ID == "13" ]; then
		echo; echo `date` "---------------------------FILTER13: Remove sites with mean DP over all individuals greater than X  -----------------------------"
		THRESHOLD=($(grep -P '^\t* *13\t* *vcftools\t* *--max-meanDP' ${CONFIG_FILE} | sed 's/\t* *13\t* *vcftools\t* *--max-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=250; fi
		#calculate the depth across loci with VCFtools
		vcftools --vcf $VCF_FILE --site-depth --out $DataName.$CutoffCode.${MODE}.Filter13
		#Now let’s take VCFtools output and cut it to only the depth scores
		cut -f3 $DataName.$CutoffCode.${MODE}.Filter13.ldepth > $DataName.$CutoffCode.${MODE}.Filter13.site.depth
		#Now let’s calculate the average depth by dividing the above file by the number of individuals
		NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $VCF_FILE )
		mawk '!/D/' $DataName.$CutoffCode.${MODE}.Filter13.site.depth | mawk -v x=$NumInd '{print $1/x}' > $DataName.$CutoffCode.${MODE}.Filter13.meandepthpersite
		cp $DataName.$CutoffCode.${MODE}.Filter13.meandepthpersite meandepthpersite
		
gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale
	set xrange [0:100] 
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
	set xrange [0:500] 
	unset label
	set title "Histogram of mean depth per site"
	set ylabel "Number of Occurrences"
	set xlabel "Mean Depth"
	binwidth=5
	bin(x,width)=width*floor(x/width) + binwidth/2.0
	set xtics 25
	plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
	pause -1
EOF
   
	gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale
	set xrange [0:1000] 
	unset label
	set title "Histogram of mean depth per site"
	set ylabel "Number of Occurrences"
	set xlabel "Mean Depth"
	binwidth=10
	bin(x,width)=width*floor(x/width) + binwidth/2.0
	set xtics 50
	plot 'meandepthpersite' using (bin($1,binwidth)):(1.0) smooth freq with boxes
	pause -1
EOF
   
		#Loci that have high mean depth are indicative of either paralogs or multicopy loci. Either way we want to remove them. Here, 
		# Remove sites with mean DP over all individuals greater than X 
		vcftools --vcf $VCF_FILE --max-meanDP ${THRESHOLD} --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter13
		mv meandepthpersite $DataName.$CutoffCode.${MODE}.Filter13.meandepthpersite

	elif [ $FILTER_ID == "14" ]; then
		echo; echo `date` "---------------------------FILTER14: If individual's genotype has DP < X, convert to missing data -----------------------------"
		THRESHOLD=($(grep -P '^\t* *14\t* *vcftools\t* *--minDP' ${CONFIG_FILE} | sed 's/\t* *14\t* *vcftools\t* *--minDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=3; fi
		vcftools --vcf $VCF_FILE --minDP ${THRESHOLD} --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter14

	elif [ $FILTER_ID == "15" ]; then
		echo; echo `date` "---------------------------FILTER15: Remove sites with maf > minor allele frequency > max-maf -----------------------------"
		THRESHOLDa=($(grep -P '^\t* *15\t* *vcftools\t* *--maf' ${CONFIG_FILE} | sed 's/\t* *15\t* *vcftools\t* *--maf\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDa}" ]; then ${THRESHOLDa}=0.005; fi
		THRESHOLDb=($(grep -P '^\t* *15\t* *vcftools\t* *--max-maf' ${CONFIG_FILE} | sed 's/\t* *15\t* *vcftools\t* *--max-maf\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDb}" ]; then ${THRESHOLDb}=0.995; fi
		# Remove sites with minor allele frequency: maf < x < max-maf
		# inspect the AF values in the vcf.  This will affect the frequency of rare variants
		vcftools --vcf $VCF_FILE --maf ${THRESHOLDa} --max-maf ${THRESHOLDb} --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter15

	elif [ $FILTER_ID == "16" ]; then
		echo; echo `date` "---------------------------FILTER16: Remove Individuals with Too Much Missing Data-----------------------------"
		THRESHOLD=($(grep -P '^\t* *16\t* *vcftools\t* *--missing-indv' ${CONFIG_FILE} | sed 's/\t* *16\t* *vcftools\t* *--missing-indv\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.5; fi
		vcftools --vcf $VCF_FILE --missing-indv
		rename out. $DataName.$CutoffCode.${MODE}.Filter16.out. out.imiss
		echo; echo "Missing Data Report, file=*out.imiss, Numbers Near 0 are Good"
		#cat $DataName.$CutoffCode.${MODE}.Filter16.out.imiss
		#graph missing data for individuals
		mawk '!/IN/' $DataName.$CutoffCode.${MODE}.Filter16.out.imiss | cut -f5 > $DataName.$CutoffCode.${MODE}.Filter16.totalmissing
		cp $DataName.$CutoffCode.${MODE}.Filter16.totalmissing totalmissing
gnuplot << \EOF 
		#filename=system("echo $DataName.$CutoffCode.${MODE}.Filter16.totalmissing")
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
		mawk '$5 > ${THRESHOLD}' $DataName.$CutoffCode.${MODE}.Filter16.out.imiss | cut -f1 > $DataName.$CutoffCode.${MODE}.Filter16.lowDP-2.indv
		##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		#list of individuals to remove
		echo `date` " Individuals with too much missing data:"
		cat $DataName.$CutoffCode.${MODE}.Filter16.lowDP-2.indv

		#remove individuals with low reads
		vcftools --vcf $VCF_FILE --remove $DataName.$CutoffCode.${MODE}.Filter16.lowDP-2.indv --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter16                                      
		mv totalmissing $DataName.$CutoffCode.$MODE.Filter16.totalmissing
		
	elif [ $FILTER_ID == "17" ]; then
		echo; echo `date` "---------------------------FILTER17: Remove sites with data missing for too many individuals in a population -----------------------------"
		THRESHOLD=($(grep -P '^\t* *17\t* *vcftools\t* *--missing-sites' ${CONFIG_FILE} | sed 's/\t* *17\t* *vcftools\t* *--missing-sites\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.5; fi
		# restrict the data to loci with a high percentage of individuals that geneotyped
		#SitesMissData: on/off switch for this filter at beginning of script
		echo "Using PopMap File: $PopMap"
		#cat $PopMap
		# get the popnames from the popmap
		popnames=`mawk '{print $2}' $PopMap | sort | uniq `
		missingDataByPop() {
			mawk -v pop=$1 '$2 == pop' $2 > $4.$3.$5.$1.keep
			vcftools --vcf $4.$3.Filter16.recode.vcf --keep $4.$3.$5.$1.keep --missing-site --out $4.$3.$5.$1
		}
		export -f missingDataByPop
		parallel --no-notice "missingDataByPop {} $PopMap $CutoffCode $DataName ${MODE}" ::: ${popnames[*]}

		# remove loci with more than X proportion of missing data
		# if you have mixed sequence lengths, this will affect if the longer regions are typed
		# UPDATE % missing threshold here	
		cat $DataName.$CutoffCode.${MODE}.*.lmiss | mawk '!/CHR/' | mawk '$6 > ${THRESHOLD}' | cut -f1,2 | sort | uniq > $DataName.$CutoffCode.${MODE}.badloci
		vcftools --vcf $VCF_FILE --exclude-positions $DataName.$CutoffCode.${MODE}.badloci --recode --recode-INFO-all --out $DataName.$CutoffCode.${MODE}.Filter17

	elif [ $FILTER_ID == "18" ]; then
		echo; echo `date` "---------------------------FILTER18: Remove sites not in HWE p<X) -----------------------------"
		THRESHOLD=($(grep -P '^\t* *18\t* *filter_hwe_by_pop_HPC' ${CONFIG_FILE} | sed 's/\t* *18\t* *filter_hwe_by_pop_HPC\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.001; fi
		# Typically, errors would have a low p-value (h setting) and would be present in many populations.
		perl $HWE_SCRIPT -v $VCF_FILE -p $PopMap -h ${THRESHOLD} -d $DataName -co $CutoffCode -o $DataName.$CutoffCode.$MODE.Filter18.HWE
		mawk '!/#/' $DataName.$CutoffCode.$MODE.Filter18.HWE.recode.vcf | wc -l
		
	elif [ $FILTER_ID == "19" ]; then
		echo; echo `date` "---------------------------FILTER19: Run rad_haplotyper to id paralogs,create haplotypes, etc -----------------------------"
		THRESHOLDa=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-d\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDa}" ]; then ${THRESHOLDa}=50; fi
		THRESHOLDb=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-mp' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-mp\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDb}" ]; then ${THRESHOLDb}=10; fi
		THRESHOLDc=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-u' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-u\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDc}" ]; then ${THRESHOLDd}=30; fi
		THRESHOLDd=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-ml' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-ml\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDd}" ]; then ${THRESHOLDd}=10; fi
		THRESHOLDe=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-h' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-h\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDe}" ]; then ${THRESHOLDe}=100; fi
		THRESHOLDf=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-z' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-z\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDf}" ]; then ${THRESHOLDf}=0.1; fi
		THRESHOLDg=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-m' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-m\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDg}" ]; then ${THRESHOLDg}=0.5; fi
		# restrict the data to loci with a high percentage of individuals that geneotyped
		###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
		#note, rad haplotyper skips complex polymorphisms by default
		#these are the Puritz Filtering Tutorial default settings
		#perl rad_haplotyper115HPC.pl -v $DataName.$CutoffCode.Filter16.HWE.recode.vcf -x $NumProc -mp 1 -u 20 -ml 4 -n -r $MAPdir/reference.$CutoffCode.fasta -bp $MAPdir -co $CutoffCode -dn $DataName  -p $MAPdir/popmap.$CutoffCode$Pops -o $DataName.$CutoffCode.Filter17.Haplotyped.vcf -g $DataName.$CutoffCode$Pops.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima  
	
		#perl rad_haplotyper116HPC.pl -v $DataName.$CutoffCode.Filter18.HWE.recode.vcf -x $NumProc -e -d 50 -mp 10 -u 30 -ml 10 -h 100 -z 0.2 -m 0.6 -r $MAPdir/reference.$CutoffCode.fasta -bp $BAM_PATH -co $CutoffCode -dn $DataName  -p $PopMap -o $DataName.$CutoffCode.Filter19.Haplotyped.vcf -g $DataName.$CutoffCode$PopMap.haps.genepop -a $DataName.$CutoffCode$Pops.haps.ima
		perl ${RADHAP_SCRIPT} -v $VCF_FILE -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -co ${CutoffCode} -dn ${DataName}  -p ${PopMap} -o ${DataName}.${CutoffCode}.${MODE}.Filter19.Haplotyped.vcf -g ${DataName}.${CutoffCode}.${MODE}.${PopMap}.haps.genepop -a ${DataName}.${CutoffCode}.${MODE}.${PopMap}.haps.ima

		#Why and how many loci were removed by rad_haplotyper
		Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
		parallel --gnu --null "grep {} $DataName.$CutoffCode.stats.out | cut -f1 > $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
		parallel --gnu --null "echo -n {}' removed, ' && wc -l $DataName.$CutoffCode.{}" ::: "${Explanations[@]}"
		echo file ${DataName}.${CutoffCode}.${MODE}.Filter19.Haplotyped.vcf has 
		mawk '!/#/' ${DataName}.${CutoffCode}.${MODE}.Filter19.Haplotyped.vcf | wc -l
		echo SNPs

	elif [ $FILTER_ID == "20" ]; then
		echo; echo `date` "---------------------------FILTER20: Select 1 Random SNP per Contig -----------------------------"
		#Script to take random SNP from every contig in a vcffile
		#Get name of VCF file
		NAME=$(echo $VCF_FILE | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 
		#Calculate number of SNPs
		Loci=(`mawk '!/#/' $VCF_FILE | wc -l `)
		#Generate list of random numbers
		seq 1 500000 | shuf | head -$Loci > nq
		#create temporary file that has a random number assigned to each SNP in first column
		cat <(mawk '/^#/' $VCF_FILE) <(paste <(mawk '!/#/' $VCF_FILE | cut -f1-5) nq <(mawk '!/#/' $VCF_FILE | cut -f7- ) )> $NAME.1RandSNP.temp
		#Use awk (mawk) to parse file and select one snp per contig (one with largest random number)
		cat $NAME.1RandSNP.temp | mawk 'BEGIN{last_loc = 0} { 
			if ($1 ~/#/) print $0;
			else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
			else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
			else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' | mawk 'NF > 0' > $NAME.randSNPperLoc.vcf
		#Remove temp file
		rm $NAME.1RandSNP.temp
		mawk '!/#/' $NAME.randSNPperLoc.vcf  | wc -l 
		mv nq ${NAME}.nq
		
	fi
}

###################################################################################################################
#help info / manual
###################################################################################################################

NAME="$(basename "$0") v1  -- a program to filter vcf files with RAD data"

SYNOPSIS="$(basename "$0") [filter settings] [input files] [output file prefix] [number of threads]"

read -d '' DESCRIPTION <<"BLOCK"
fltrVCF is a tool to filter VCF files created by dDocentHPC. BfltrVCF is designed to be
        run in stages, e.g., cleaning then polishing. Cleaning attemps to keep biological variation
        and remove erroneous and poorly supported variation.  Polishing removes undesireable
        variation that is biological, e.g., loci or samples with too much missing data,
        indels, SNPs with >2 alleles, etc. In either stage, you can run any of the filters in any
        order you choose.

        The default filter settings files are either config.fltr.clean or config.fltr.polish and
        are not interchangeable. To change the settings, copy and rename the config files so that
        you don't lose the defaults.

        fltrVCF is parallelized where possible, but only runs on one node or computer. MPI is not
        supported.

        fltrVCF requires minor modification to work with dDocent output.  To do so, remove
        ".$CutoffCode" "$CutoffCode." and "$CutoffCode" in order from the script). Both
        filter_hwe_by_pop_HPC.pl and rad_haplotyperHPC116.pl are tested with fltrVCF and work. Use of
        other versions is possible, and will be neccessary if filtering data created by dDocent rather
        than dDocentHPC, but is not supported.
BLOCK

read -d '' OPTIONS <<"BLOCK"
[filter settings]
                -m <arg>        filtering stage, a label that will be included with output files.
                                 This allows data to be processed in stages in the same directory
                                 without overwriting. For example, -m can be set to "clean" for the
                                 first stage of filtering, where errors are filtered and all variation
                                 perceived to be biological is retained. In the same directory, a
                                 second phase of filtering can be run without overwriting previous
                                 output by setting -m to be "polish". In this second phase, variation
                                 perceived to be biological, but undesireable is filtered to obtain
                                 a desireable data set. This value can be set to any string without
                                 white space. [clean]
                -f <arg>        if set, controls filters to be run, in order. Argument should be 2
                                 digit numbers separated by spaces. -f "01 04 02"  or  -f 01\ 04\ 02
                                 will specify that filters 01, 04, and 02 will be run in succession.
                                 Filters are described in the config files. If -f is not set, the
                                 config file is used to determine the filters and order. If -f is
                                 set, it will override the config file. []
                -s <arg>        file with filter settings [config.fltr.clean.ind]
        [input files]
                -c <arg>        cutoff values used for reference genome [3.3]
                -b <arg>        path to mapping directory with *.bam [../mapping]
                -v <arg>        vcf file to be filtered [${b}/TotalRawSNPs.${c}.vcf]
                -g <arg>        reference genome fasta file [${b}/reference.${c}.fasta]
                -p <arg>        popmap file to use for defining population affiliation
                                 [${b}/popmap.${c}]
                -w <arg>        filter_hwe perl script [filter_hwe_by_pop_HPC.pl]
                -r <arg>        rad_haplotyper perl script [rad_haplotyperHPC116.pl]

        [output file prefix]
                -o <arg>        optional, all output files will be prefixed with this argument []

        [number of threads]
                -t <arg>        number of threads available for parallel processing [1]

EXAMPLES
        The following two commands are the same, the first takes advantage of the defaults,
        the second does not.

                fltrVCF.bash -c 25.10 -o ProjectX.A -t 40

                fltrVCF.bash -m clean -I -c 25.10 -m ../mapping -v ../mapping/TotalRawSNPs.3.6.vcf
                        -p ../mapping/popmap.25.10 -s config.fltr.clean -w filter_hwe_by_pop.pl
                        -r rad_haplotyperHPC116.pl -o ProjectX.A -t 40

BLOCK


###################################################################################################################
#Read in command line arguments
###################################################################################################################
while getopts ":m:f:c:b:v:g:p:s:w:r:o:t:h" opt; do
  case $opt in
    h)
        echo ""
        echo "NAME" >&2
        echo -e '\t'"$NAME" >&2
        echo ""
        echo "SYNOPSIS" >&2
        echo -e '\t'"$SYNOPSIS"
        echo ""
        echo "DESCRIPTION"
        echo -e '\t'"$DESCRIPTION" >&2
        echo ""
        echo "OPTIONS" >&2
        echo -e '\t'"$OPTIONS" >&2
        echo ""
      exit
      ;;
    m)
        if [[ ${OPTARG} == "clean" || ${OPTARG} == "polish" ]]; then
                echo "Mode:                 $OPTARG"
                MODE=$OPTARG
        else
                echo "ERROR :-(                 Invalid -m option: $OPTARG" >&2
                exit
        fi
        ;;
	\?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    f)
        echo "Filter order:                 $OPTARG"
        FILTERS=($OPTARG)
        ;;
	\?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    c)
        echo "" >&2
        echo "Cutoffs:                  $OPTARG"
        CutoffCode=$OPTARG >&2
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    v)
        echo "VCF File:                 $OPTARG"
        VCF=$OPTARG
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    g)
        echo "REF File:                 $OPTARG"
        REF_FILE=$OPTARG
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
	b)
        echo "Path to BAM files:        $OPTARG"
        BAM_PATH=$OPTARG
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    p)
        echo "PopMap File:              $OPTARG"
        PopMap=$OPTARG
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    s)
        echo "Settings File:            $OPTARG"
        CONFIG_FILE=$OPTARG
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    w)
        echo "HWE Script:               $OPTARG"
        HWE_SCRIPT=$(${OPTARG})
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    r)
        echo "Rad_Haplotyper script:    $OPTARG"
        RADHAP_SCRIPT=$(${OPTARG})
        ;;
    \?)
        echo "ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    o)
        echo "Output file prefix:       $OPTARG"
        DataName=$OPTARG
        ;;
    \?)
      echo "ERROR :-/                   Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "ERROR :-/                   Option -$OPTARG requires an argument." >&2
      ;;
    t)
      echo "Number of threads:          $OPTARG"
      ;;
    \?)
      echo "ERROR :-/                   Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "ERROR :-/                   Option -$OPTARG requires an argument." >&2
      ;;
  esac
done

echo "" >&2


###################################################################################################################
#Set defaults if user input and config file info is absent
###################################################################################################################

if [ -z ${MODE+x} ]; then MODE="clean" && echo "Filtering phase is set to default '${MODE}'"; else echo "Filtering phase is set to '${MODE}'"; fi

if [ -z ${CONFIG_FILE+x} ]; then 
	CONFIG_FILE=$(ls config.fltr.clean.ind)
	echo "Settings are being loaded from: '${CONFIG_FILE}' by default"
else 
	CONFIG_FILE=$(ls ${CONFIG_FILE}) 
	echo "SETTINGS are being loaded from file: '${CONFIG_FILE}'"
fi
if [ $CONFIG_FILE == "" ]; then 
	echo "ERROR :-< 	configuration file is missing" >&2
	exit
fi

if [ -z ${FILTERS+x} ]; then 
	FILTERS=($(grep -P '^\t* *fltrVCF\t* *-f\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	echo "Filters are set to '${FILTERS[@]}'"
else 
	echo "FILTERS are set to '${FILTERS[@]}'"
fi
if [ ${FILTERS[0]} == "" ]; then 
	echo "ERROR :-< 	filter settings (-f) are missing.  Please modify -f argument in config file or at the command line" >&2
	exit
fi

if [ -z ${CutoffCode+x} ]; then 
	CutoffCode=($(grep -P '^\t* *fltrVCF\t* *-c\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $CutoffCode == "" ]; then
		CutoffCode="3.3" 
	fi
fi
echo "CutoffCode is set to '${CutoffCode}'"

if [ -z ${BAM_PATH+x} ]; then 
	BAM_PATH=($(grep -P '^\t* *fltrVCF\t* *-b\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $BAM_PATH == "" ]; then
		BAM_PATH="../mapping" 
	fi
fi
echo "BAM_PATH is set to '${BAM_PATH}'"

if [ -z ${VCF_FILE+x} ]; then 
	VCF_FILE=($(grep -P '^\t* *fltrVCF\t* *-v\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $VCF_FILE == "" ]; then
		VCF_FILE=${BAM_PATH}"/TotalRawSNPs."${CutoffCode}".vcf"
	fi
fi
echo "VCF_FILE is set to '${VCF_FILE}'"

if [ -z ${REF_FILE+x} ]; then 
	REF_FILE=($(grep -P '^\t* *fltrVCF\t* *-g\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $REF_FILE == "" ]; then
		REF_FILE=${BAM_PATH}"/reference."${CutoffCode}".fasta" 
	fi
fi	
echo "Reference genome is set to '${REF_FILE}'"

if [ -z ${PopMap+x} ]; then 
	PopMap=($(grep -P '^\t* *fltrVCF\t* *-p\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $PopMap == "" ]; then
		PopMap=${BAM_PATH}"/popmap."${CutoffCode} 
	fi
fi
echo "PopMap is set to '${PopMap}'"

if [ -z ${HWE_SCRIPT+x} ]; then 
	HWE_SCRIPT=($(grep -P '^\t* *fltrVCF\t* *-w\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $HWE_SCRIPT == "" ]; then
		HWE_SCRIPT="filter_hwe_by_pop_HPC.pl" 
	fi
fi
echo "HWE_SCRIPT is set to '${HWE_SCRIPT}'"

if [ -z ${RADHAP_SCRIPT+x} ]; then 
	RADHAP_SCRIPT=($(grep -P '^\t* *fltrVCF\t* *-r\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $RADHAP_SCRIPT == "" ]; then
		RADHAP_SCRIPT="rad_haplotyperHPC116.pl" 
	fi
fi
echo "RADHAP_SCRIPT is set to '${RADHAP_SCRIPT}'"

if [ -z ${DataName+x} ]; then 
	DataName=($(grep -P '^\t* *fltrVCF\t* *-o\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $DataName == "" ]; then
		echo "No output file prefix is specified.'" 
	else
		echo "Output file prefix is set to '${DataName}'"
	fi
else
	echo "Output file prefix is set to '${DataName}'"
fi


if [ -z ${NumProc+x} ]; then 
	NumProc=($(grep -P '^\tfltrVCF -t ' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $NumProc == "" ]; then
		NumProc=1 
	fi
fi
echo "The number of threads is set to '${NumProc}'"


###################################################################################################################
#Get filter thresholds from config file
###################################################################################################################



###################################################################################################################
#Run script
###################################################################################################################

MAIN FILTERS $MODE $CutoffCode $BAM_PATH $VCF_FILE $REF_FILE $PopMap $CONFIG_FILE $HWE_SCRIPT $RADHAP_SCRIPT $DataName $NumProc
 