#!/bin/bash
VERSION=3.2
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
	# echo `date` " reading variables into MAIN"
	name=$1[@]
	FILTERS=("${!name}")
	CutoffCode=$2
	BAM_PATH=$3
	VCF_FILE=$4
	REF_FILE=$5
	PopMap=$6
	CONFIG_FILE=$7
	HWE_SCRIPT=$8
	RADHAP_SCRIPT=$9
	DataName=${10}
	NumProc=${11}
	PARALLEL=${12}
	
	# echo "	${FILTERS[0]}"
	# echo "	$CutoffCode"
	# echo $BAM_PATH
	# echo $VCF_FILE
	# echo $REF_FILE
	# echo $PopMap
	# echo $CONFIG_FILE
	# echo $HWE_SCRIPT
	# echo $RADHAP_SCRIPT
	# echo $DataName
	# echo $NumProc
	# echo $PARALLEL
	
	for i in ${FILTERS[@]}; do
		FILTER $i $CutoffCode $BAM_PATH $VCF_FILE $REF_FILE $PopMap $CONFIG_FILE $HWE_SCRIPT $RADHAP_SCRIPT $DataName $NumProc $PARALLEL
		echo $i
		if [ $PARALLEL == "FALSE" ]; then
			VCF_FILE=$(ls -t ${DataName}*${i}*vcf | head -n 1)
		else
			VCF_FILE=$(ls -t ${DataName}*${i}*vcf.gz | head -n 1)
		fi
	done

}


###################################################################################################################
#Define filters 
###################################################################################################################
function FILTER(){

#	echo `date` " reading in variables to FILTER"
	FILTER_ID=$1
	CutoffCode=$2
	BAM_PATH=$3
	VCF_FILE=$4
	REF_FILE=$5
	PopMap=$6
	CONFIG_FILE=$7
	HWE_SCRIPT=$8
	RADHAP_SCRIPT=$9
	DataName=${10}
	NumProc=${11}
	PARALLEL=${12}

	# echo $FILTER_ID
	# echo $CutoffCode
	# echo $BAM_PATH
	# echo $VCF_FILE
	# echo $REF_FILE
	# echo $PopMap
	# echo $CONFIG_FILE
	# echo $HWE_SCRIPT
	# echo $RADHAP_SCRIPT
	# echo $DataName
	# echo $NumProc
	# echo $PARALLEL
	
	
	if [ $FILTER_ID == "01" ]; then
		echo; echo `date` "---------------------------FILTER01: Remove sites with Y < alleles < X -----------------------------"
		#get settings from config file
		THRESHOLD=($(grep -P '^\t* *01\t* *vcftools\t* *--min-alleles' ${CONFIG_FILE} | sed 's/\t* *01\t* *vcftools\t* *--min-alleles\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "$THRESHOLD" ]; then THRESHOLD=2; fi
		THRESHOLDb=($(grep -P '^\t* *01\t* *vcftools\t* *--max-alleles' ${CONFIG_FILE} | sed 's/\t* *01\t* *vcftools\t* *--max-alleles\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "$THRESHOLDb" ]; then THRESHOLDb=2; fi
		Filter="\"--min-alleles $THRESHOLD --max-alleles $THRESHOLDb --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		#call function to filter vcf
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc

	elif [ $FILTER_ID == "02" ]; then
		echo; echo `date` "---------------------------FILTER02: Remove Sites with Indels -----------------------------"
		Filter="\"--remove-indels --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "03" ]; then
		echo; echo `date` "---------------------------FILTER03: Remove Sites with QUAL < minQ -----------------------------"
		THRESHOLD=($(grep -P '^\t* *03\t* *vcftools\t* *--minQ' ${CONFIG_FILE} | sed 's/\t* *03\t* *vcftools\t* *--minQ\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=30; fi
		Filter="\"--minQ ${THRESHOLD} --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "04" ]; then
		echo; echo `date` "---------------------------FILTER04: Remove Sites With Mean Depth of Coverage < min-meanDP -----------------------------"
		THRESHOLD=($(grep -P '^\t* *04\t* *vcftools\t* *--min-meanDP' ${CONFIG_FILE} | sed 's/\t* *04\t* *vcftools\t* *--min-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=2; fi
		Filter="\"--min-meanDP ${THRESHOLD} --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "05" ]; then
		echo; echo `date` "---------------------------FILTER05: Remove sites called in <X proportion of individuals -----------------------------"
		THRESHOLD=($(grep -P '^\t* *05\t* *vcftools\t* *--max-missing' ${CONFIG_FILE} | sed 's/\t* *05\t* *vcftools\t* *--max-missing\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.5; fi
		Filter="\"--max-missing ${THRESHOLD} --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "06" ]; then
		echo; echo `date` "---------------------------FILTER06: Remove sites with Average Allele Balance deviating too far from 0.5 while keeping those with AB=0  -----------------------------"
		THRESHOLD=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *min' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.375; fi
		THRESHOLDb=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *max' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *max\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDb}" ]; then ${THRESHOLDb}=0.625; fi
		Filter="\"AB > $THRESHOLD & AB < $THRESHOLDb | AB = 0\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "TRUE" #the last option "TRUE" is for vcffixup

#UNDER CONSTRUCTION V
	elif [ $FILTER_ID == "96" ]; then
		echo; echo `date` "---------------------------FILTER06: Remove contigs with Average Allele Balance deviating too far from 0.5 while keeping those with AB=0  -----------------------------"
		THRESHOLD=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *min' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.375; fi
		THRESHOLDb=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *max' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *max\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDb}" ]; then ${THRESHOLDb}=0.625; fi
		Filter="\"AB > $THRESHOLD & AB < $THRESHOLDb | AB = 0\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		#get the AB data
		grep '^dDocent_Contig' $VCF_FILE | cut -f1,8 | sed 's/;/\t/g' | cut -f1,2 | sed 's/AB=//' > ${VCF_OUT%.*}.AB.txt
		#calculate mean AB
		awk '{print $1}' AB.txt | uniq | uniq.contigs.txt
		awk '$1 == "dDocent_Contig_1"' AB.txt | awk '$2 != 0' | awk '{total += $2; n++ } END {if (n > 0) print "dDocent_Contig_1 " total / n} '
		#do it in parallel
		#cat uniq.contigs.txt | parallel --no-notice -k -j 
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "TRUE" #the last option "TRUE" is for vcffixup
#UNDER CONSTRUCTION ^
	elif [ $FILTER_ID == "07" ]; then
		echo; echo `date` "---------------------------FILTER07: Remove sites with Alternate Allele Count <=X -----------------------------"
		THRESHOLD=($(grep -P '^\t* *07\t* *vcffilter\t* *AC\t* *min' ${CONFIG_FILE} | sed 's/\t* *07\t* *vcffilter\t* *AC\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=1; fi
		Filter="\"AC > $THRESHOLD & AN - AC > $THRESHOLD\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [ $FILTER_ID == "08" ]; then
		echo; echo `date` "---------------------------FILTER08: Remove sites covered by both F and R reads-----------------------------"
		# This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)
		THRESHOLD=($(grep -P '^\t* *08\t* *vcffilter\t* *SAF\/SAR\t* *min' ${CONFIG_FILE} | sed 's/\t* *08\t* *vcffilter\t* *SAF\/SAR\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=1; fi
		Filter="\"SAF / SAR > ${THRESHOLD} & SRF / SRR > ${THRESHOLD} | SAR / SAF > ${THRESHOLD} & SRR / SRF > ${THRESHOLD}\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [ $FILTER_ID == "09" ]; then
		echo; echo `date` "---------------------------FILTER09: Remove sites with low/high ratio of mean mapping quality of alt to ref -----------------------------"
		# The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
		# loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
		THRESHOLD=($(grep -P '^\t* *09\t* *vcffilter\t* *MQM\/MQMR\t* *min' ${CONFIG_FILE} | sed 's/\t* *09\t* *vcffilter\t* *MQM\/MQMR\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.1; fi
		calc(){ awk "BEGIN { print "$*" }"; }
		THRESHOLDmin=$(calc 1-${THRESHOLD})
		#THRESHOLDmin=$(echo 1-${THRESHOLD} | R --vanilla --quiet | sed -n '2s/.* //p')
		THRESHOLDmax=$(calc 1/${THRESHOLDmin})
		#THRESHOLDmax=$(echo 1/${THRESHOLDmin} | R --vanilla --quiet | sed -n '2s/.* //p')
		Filter="\"MQM / MQMR > 1 - ${THRESHOLDmin} & MQM / MQMR < ${THRESHOLDmax}\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [ $FILTER_ID == "10" ]; then
		echo; echo `date` "---------------------------FILTER10: Remove sites with one allele only supported by reads not properly paired -----------------------------"
		# another filter that can be applied is whether or not their is a discrepancy in the properly paired status  for reads supporting reference or alternate alleles.
		# Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
		# the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
		Filter="\"PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"
		
	elif [ $FILTER_ID == "11" ]; then
		echo; echo `date` "---------------------------FILTER11: Remove sites with Quality/DP ratio < 0.25 -----------------------------"
		# There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.
		# first, by removing any locus that has a quality score below 1/4 of the read depth.
		THRESHOLD=($(grep -P '^\t* *11\t* *vcffilter\t* *QUAL\/DP\t* *min' ${CONFIG_FILE} | sed 's/\t* *11\t* *vcffilter\t* *QUAL\/DP\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.25; fi
		Filter="\"QUAL / DP > ${THRESHOLD}\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [ $FILTER_ID == "12" ]; then
		echo; echo `date` "---------------------------FILTER12: Remove sites with Quality > 2xDP -----------------------------"
		# second, is a multistep process
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		Filter="\"--exclude-positions $VCF_OUT.lowQDloci --recode --recode-INFO-all\""
		# create a list of the depth of each locus
		cut -f8 $VCF_FILE | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $VCF_OUT.DEPTH
		#Then make a list of quality scores
		mawk '!/#/' $VCF_FILE | cut -f1,2,6 > $VCF_OUT.loci.qual
		#Then calculate mean depth
		MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $VCF_OUT.DEPTH )
		#echo $MeanDepth
		MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
		#echo $MeanDepth2
		#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
		paste $VCF_OUT.loci.qual $VCF_OUT.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $VCF_OUT.lowQDloci
		#cat $VCF_OUT.lowQDloci
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "13" ]; then
		echo; echo `date` "---------------------------FILTER13: Remove sites with mean DP greater than X  -----------------------------"
		THRESHOLD=($(grep -P '^\t* *13\t* *vcftools\t* *--max-meanDP' ${CONFIG_FILE} | sed 's/\t* *13\t* *vcftools\t* *--max-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=250; fi
		Filter="--site-depth"
		Filter2="\"--max-meanDP $THRESHOLD --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		#calculate the depth across loci with VCFtools
		if [ $PARALLEL == "FALSE" ]; then 
			vcftools --vcf ${VCF_FILE} $Filter --out $VCF_OUT 2> /dev/null
			#Now let’s take VCFtools output and cut it to only the depth scores
			cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			#Now let’s calculate the average depth by dividing the above file by the number of individuals
			NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $VCF_FILE )
		else
			vcftools --gzvcf ${VCF_FILE} $Filter --out $VCF_OUT 2> /dev/null
			#Now let’s take VCFtools output and cut it to only the depth scores
			cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			#Now let’s calculate the average depth by dividing the above file by the number of individuals
			NumInd=$(gunzip -c $VCF_FILE | awk '{if ($1 == "#CHROM"){print NF-9; exit}}' )
		fi
		mawk '!/D/' $VCF_OUT.site.depth | mawk -v x=$NumInd '{print $1/x}' > $VCF_OUT.meandepthpersite
		cp $VCF_OUT.meandepthpersite meandepthpersite
		
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
   
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter2}" $VCF_OUT $DataName $CutoffCode $NumProc 
		rm meandepthpersite

	elif [ $FILTER_ID == "14" ]; then
		echo; echo `date` "---------------------------FILTER14: If individual's genotype has DP < X, convert to missing data -----------------------------"
		THRESHOLD=($(grep -P '^\t* *14\t* *vcftools\t* *--minDP' ${CONFIG_FILE} | sed 's/\t* *14\t* *vcftools\t* *--minDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=3; fi
		Filter="\"--minDP ${THRESHOLD} --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "15" ]; then
		echo; echo `date` "---------------------------FILTER15: Remove sites with maf > minor allele frequency > max-maf -----------------------------"
		THRESHOLDa=($(grep -P '^\t* *15\t* *vcftools\t* *--maf' ${CONFIG_FILE} | sed 's/\t* *15\t* *vcftools\t* *--maf\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDa}" ]; then ${THRESHOLDa}=0.005; fi
		THRESHOLDb=($(grep -P '^\t* *15\t* *vcftools\t* *--max-maf' ${CONFIG_FILE} | sed 's/\t* *15\t* *vcftools\t* *--max-maf\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDb}" ]; then ${THRESHOLDb}=0.995; fi
		# Remove sites with minor allele frequency: maf < x < max-maf
		# inspect the AF values in the vcf.  This will affect the frequency of rare variants
		Filter="\"--maf ${THRESHOLDa} --max-maf ${THRESHOLDb} --recode --recode-INFO-all\""
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "16" ]; then
		echo; echo `date` "---------------------------FILTER16: Remove Individuals with Too Much Missing Data-----------------------------"
		THRESHOLD=($(grep -P '^\t* *16\t* *vcftools\t* *--missing-indv' ${CONFIG_FILE} | sed 's/\t* *16\t* *vcftools\t* *--missing-indv\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.5; fi
		Filter="--missing-indv"
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		if [ $PARALLEL == "FALSE" ]; then 
			vcftools --vcf ${VCF_FILE} $Filter --out $VCF_OUT 2> /dev/null
		else
			vcftools --gzvcf ${VCF_FILE} $Filter --out $VCF_OUT 2> /dev/null
		fi
		echo; echo "Missing Data Report, file=*out.imiss, Numbers Near 0 are Good"
		#cat $VCF_OUT.imiss
		#graph missing data for individuals
		mawk '!/IN/' $VCF_OUT.imiss | cut -f5 > $VCF_OUT.totalmissing
		cp $VCF_OUT.totalmissing totalmissing
gnuplot << \EOF 
		#filename=system("echo $VCF_OUT.totalmissing")
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
		mawk -v x="$THRESHOLD" '$5 > x' $VCF_OUT.imiss | cut -f1 > $VCF_OUT.lowDP-2.indv
		##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Filter2="\"--remove $VCF_OUT.lowDP-2.indv --recode --recode-INFO-all\""
		#list of individuals to remove
		echo `date` " Individuals with too much missing data:"
		cat $VCF_OUT.lowDP-2.indv
		#remove individuals with low reads
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter2}" $VCF_OUT $DataName $CutoffCode $NumProc 
		rm totalmissing
		
	elif [ $FILTER_ID == "17" ]; then
		echo; echo `date` "---------------------------FILTER17: Remove sites with data missing for too many individuals in a population -----------------------------"
		THRESHOLD=($(grep -P '^\t* *17\t* *vcftools\t* *--missing-sites' ${CONFIG_FILE} | sed 's/\t* *17\t* *vcftools\t* *--missing-sites\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.5; fi
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		# restrict the data to loci with a high percentage of individuals that geneotyped
		#SitesMissData: on/off switch for this filter at beginning of script
		echo "Using PopMap File: $PopMap"
		#cat $PopMap
		# get the popnames from the popmap
		popnames=`mawk '{print $2}' $PopMap | sort | uniq `
		missingDataByPop() {
			PARALLEL=$3
			VCF_FILE=$4
			VCF_OUT=$5
			mawk -v pop=$1 '$2 == pop' $2 > $VCF_OUT.$1.keep
			if [ $PARALLEL == "FALSE" ]; then
				vcftools --vcf $VCF_FILE --keep $VCF_OUT.$1.keep --missing-site --out $VCF_OUT.$1 2> /dev/null
			else
				vcftools --gzvcf $VCF_FILE --keep $VCF_OUT.$1.keep --missing-site --out $VCF_OUT.$1 2> /dev/null
			fi
		}
		export -f missingDataByPop
		parallel --no-notice -k -j $NumProc "missingDataByPop {} $PopMap $PARALLEL $VCF_FILE $VCF_OUT 2> /dev/null" ::: ${popnames[*]}

		# remove loci with more than X proportion of missing data
		# if you have mixed sequence lengths, this will affect if the longer regions are typed
		# UPDATE % missing threshold here	
		cat $VCF_OUT.*.lmiss | mawk '!/CHR/' | mawk -v x=$THRESHOLD '$6 > x' | cut -f1,2 | sort | uniq > $VCF_OUT.badloci
		Filter2="\"--exclude-positions $VCF_OUT.badloci --recode --recode-INFO-all\""
		FILTER_VCFTOOLS $PARALLEL $VCF_FILE "${Filter2}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [ $FILTER_ID == "18" ]; then
		echo; echo `date` "---------------------------FILTER18: Remove sites not in HWE p<X) -----------------------------"
		THRESHOLD=($(grep -P '^\t* *18\t* *filter_hwe_by_pop_HPC' ${CONFIG_FILE} | sed 's/\t* *18\t* *filter_hwe_by_pop_HPC\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLD}" ]; then ${THRESHOLD}=0.001; fi
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		if [ $PARALLEL == "TRUE" ]; then 
			gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			perl $HWE_SCRIPT -v ${VCF_FILE%.*} -p $PopMap -h ${THRESHOLD} -d $DataName -co $CutoffCode -o $VCF_OUT.HWE
			mawk '!/#/' $VCF_OUT.HWE.recode.vcf | wc -l
			bgzip -@ $NumProc -c $VCF_OUT.HWE.recode.vcf > $VCF_OUT.HWE.recode.vcf.gz
			tabix -p vcf $VCF_OUT.HWE.recode.vcf.gz
		else
			# Typically, errors would have a low p-value (h setting) and would be present in many populations.
			perl $HWE_SCRIPT -v $VCF_FILE -p $PopMap -h ${THRESHOLD} -d $DataName -co $CutoffCode -o $VCF_OUT.HWE
			mawk '!/#/' $VCF_OUT.HWE.recode.vcf | wc -l
		fi
		
		
	elif [ $FILTER_ID == "19" ]; then
		echo; echo `date` "---------------------------FILTER19: Run rad_haplotyper to id paralogs,create haplotypes, etc -----------------------------"
		THRESHOLDa=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-d\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-d\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDa}" ]; then ${THRESHOLDa}=50; fi
		THRESHOLDb=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-mp\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-mp\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDb}" ]; then ${THRESHOLDb}=10; fi
		THRESHOLDc=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-u\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-u\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDc}" ]; then ${THRESHOLDd}=30; fi
		THRESHOLDd=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-ml\t* *' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-ml\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDd}" ]; then ${THRESHOLDd}=10; fi
		THRESHOLDe=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-h\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-h\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDe}" ]; then ${THRESHOLDe}=100; fi
		THRESHOLDf=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-z\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-z\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDf}" ]; then ${THRESHOLDf}=0.1; fi
		THRESHOLDg=($(grep -P '^\t* *19\t* *rad_haplotyperHPC116\t* *-m\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyperHPC116\t* *-m\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [ -z "${THRESHOLDg}" ]; then ${THRESHOLDg}=0.5; fi
		VCF_OUT=$DataName.$CutoffCode
		echo "Read Sampling Depth:					$THRESHOLDa"
		echo "Max Paralogous Individuals: 				$THRESHOLDb"
		echo "Max Num SNPs per Contig:				$THRESHOLDc"
		echo "Max Low Cov Individuals:				$THRESHOLDd"
		echo "Max Haplotypes in XS of NumSNPs:			$THRESHOLDe"
		echo "Num or Prop of Reads With Rare Haplotypes to Ignore:	$THRESHOLDf"
		echo "Min Prop of Individuals with Haplotypes to Keep Locus:	$THRESHOLDg"
		
		# restrict the data to loci with a high percentage of individuals that geneotyped
		###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
		#note, rad haplotyper skips complex polymorphisms by default
		#runs ceb version of rad haplotyper
		# if [ $PARALLEL == "TRUE" ]; then 
			# gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			# perl $RADHAP_SCRIPT -v ${VCF_FILE%.*} -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -co ${CutoffCode} -dn ${DataName}  -p ${PopMap} -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf -g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			# mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			# bgzip -@ $NumProc -c $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf > $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
			# tabix -p vcf $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
		# else
			# perl $RADHAP_SCRIPT -v $VCF_FILE -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -co ${CutoffCode} -dn ${DataName}  -p ${PopMap} -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf -g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			# mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
		# fi

		#runs official rad haplotyper
		sed "s/\t/\.$CutoffCode\t/g" $PopMap > $PopMap.$CutoffCode
		
		if [ $PARALLEL == "TRUE" ]; then 
			gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			#modify individual names in vcf to include cutoffcode
			LineNum=$(awk '/^#CHROM/{print NR; exit}' ${VCF_FILE%.*})
			FirstCols=$(grep '^#CHROM' ${VCF_FILE%.*} | cut -f1-9)
			NewIndNames=$(grep '^#CHROM' ${VCF_FILE%.*} | cut -f10- | sed "s/\t/\.$CutoffCode\t/g" | sed "s/$/\.$CutoffCode/g")
			NewLine=$(echo $FirstCols $NewIndNames)
			sed "${LineNum}s/.*/$(echo $NewLine | sed "s/ /\t/g")/" ${VCF_FILE%.*} > ${VCF_FILE%.*}.1.vcf

			perl $RADHAP_SCRIPT -v ${VCF_FILE%.*}.1.vcf -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}.$CutoffCode -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf #-g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			
			
			
			mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			bgzip -@ $NumProc -c $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf > $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
			tabix -p vcf $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
			rm ${VCF_FILE%.*}.1.vcf
		else
			#modify individual names in vcf to include cutoffcode
			LineNum=$(awk '/^#CHROM/{print NR; exit}' $VCF_FILE)
			FirstCols=$(grep '^#CHROM' $VCF_FILE | cut -f1-9)
			NewIndNames=$(grep '^#CHROM' $VCF_FILE | cut -f10- | sed "s/\t/\.$CutoffCode\t/g" | sed "s/$/\.$CutoffCode/g")
			NewLine=$(echo $FirstCols $NewIndNames)
			sed "${LineNum}s/.*/$(echo $NewLine | sed "s/ /\t/g")/" $VCF_FILE > $VCF_FILE.1.vcf
			
			perl $RADHAP_SCRIPT -v $VCF_FILE.1.vcf -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}.$CutoffCode -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf #-g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			rm $VCF_FILE.1.vcf
		fi
		
		rm $PopMap.$CutoffCode

		#rename outputfiles
		if [ -f "fail.log" ]; then mv fail.log $VCF_OUT.Fltr$FILTER_ID.fail.log; fi
		if [ -f "haplo_dump.out" ]; then mv haplo_dump.out $VCF_OUT.Fltr$FILTER_ID.haplo_dump.out; fi
		if [ -f "hap_log.out" ]; then mv hap_log.out $VCF_OUT.Fltr$FILTER_ID.hap_log.out; fi
		if [ -f "stats.out" ]; then mv stats.out $VCF_OUT.Fltr$FILTER_ID.stats.out; fi
		if [ -f "ind_stats.out" ]; then mv ind_stats.out $VCF_OUT.Fltr$FILTER_ID.ind_stats.out; fi
		if [ -f "snp_dump.out" ]; then mv snp_dump.out $VCF_OUT.Fltr$FILTER_ID.snp_dump.out; fi
		if [ -f "allele_dump.out" ]; then mv allele_dump.out $VCF_OUT.Fltr$FILTER_ID.allele_dump.out; fi
		if [ -f "temp.vcf" ]; then mv temp.vcf $VCF_OUT.Fltr$FILTER_ID.temp.vcf; fi
		
		#Why and how many loci were removed by rad_haplotyper
		Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
		parallel --gnu --null "grep {} $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f1 > $VCF_OUT.Fltr$FILTER_ID.{}" ::: "${Explanations[@]}"
		parallel --gnu --null "echo -n {}' removed, ' && wc -l $VCF_OUT.Fltr$FILTER_ID.{}" ::: "${Explanations[@]}"
		echo file $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf has 
		mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
		echo SNPs

	elif [ $FILTER_ID == "20" ]; then
		echo; echo `date` "---------------------------FILTER20: Select 1 Random SNP per Contig -----------------------------"
		#Script to take random SNP from every contig in a vcffile
		#Get name of VCF file
		VCF_OUT=$DataName.$CutoffCode.Fltr$FILTER_ID
		if [ $PARALLEL == "TRUE" ]; then 
			gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			VCF_FILE=${VCF_FILE%.*}
		fi

		NAME=$(echo $VCF_FILE | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 
		#Calculate number of SNPs
		Loci=(`mawk '!/#/' $VCF_FILE | wc -l `)
		#Generate list of random numbers
		seq 1 500000 | shuf | head -$Loci > $VCF_OUT.nq
		#create temporary file that has a random number assigned to each SNP in first column
		cat <(mawk '/^#/' $VCF_FILE) <(paste <(mawk '!/#/' $VCF_FILE | cut -f1-5) $VCF_OUT.nq <(mawk '!/#/' $VCF_FILE | cut -f7- ) )> $VCF_OUT.1RandSNP.temp
		#Use awk (mawk) to parse file and select one snp per contig (one with largest random number)
		cat $VCF_OUT.1RandSNP.temp | mawk 'BEGIN{last_loc = 0} { 
			if ($1 ~/#/) print $0;
			else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$6}
			else if ($1 == last_loc && $6 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$6}
			else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$6}
		} END{print last_contig}' | mawk 'NF > 0' > $VCF_OUT.randSNPperLoc.vcf
		#Remove temp file
		rm $VCF_OUT.1RandSNP.temp
		mawk '!/#/' $VCF_OUT.randSNPperLoc.vcf  | wc -l 
		rm $VCF_OUT.nq
		if [ $PARALLEL == "TRUE" ]; then 
			bgzip -@ $NumProc -c $VCF_OUT.randSNPperLoc.vcf > $VCF_OUT.randSNPperLoc.vcf.gz
			tabix -p vcf $VCF_OUT.randSNPperLoc.vcf.gz
		fi
	fi
}


###################################################################################################################
#specific filter functions
###################################################################################################################

function FILTER_VCFTOOLS(){
	PARALLEL=$1
	VCF_FILE=$2
	Filter=$(echo "$3")
	VCF_OUT=$4
	DataName=$5
	CutoffCode=$6
	NumProc=$7
	
	if [ $PARALLEL == "FALSE" ]; then 
		vcftools --vcf $VCF_FILE $Filter --out $VCF_OUT 2> /dev/null
		echo -n "	Sites remaining:	" && mawk '!/#/' $VCF_OUT.recode.vcf | wc -l
		echo -n "	Contigs remaining:	" && mawk '!/#/' $VCF_OUT.recode.vcf | cut -f1 | uniq | wc -l
	else
		tabix -H $VCF_FILE > $VCF_OUT.header.vcf
		NumHeaderLines=$(tabix -H $VCF_FILE | wc -l)
		NumHeaderLines=$(($NumHeaderLines+1))
		ls $DataName.$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcftools --vcf - \"$Filter\" --stdout 2> /dev/null | tail -n +$NumHeaderLines" 2> /dev/null | cat $VCF_OUT.header.vcf - | bgzip -@ $NumProc -c > $VCF_OUT.recode.vcf.gz
		tabix -f -p vcf $VCF_OUT.recode.vcf.gz
		echo -n "	Sites remaining:	" && ls $DataName.*.$CutoffCode.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.recode.vcf.gz | wc -l " | awk -F: '{a+=$1} END{print a}'
		echo -n "	Contigs remaining:	" && ls $DataName.*.$CutoffCode.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.recode.vcf.gz | cut -f1 | uniq | wc -l " | awk -F: '{a+=$1} END{print a}'
		rm $VCF_OUT.header.vcf
	fi
	
}

function FILTER_VCFFILTER(){
	PARALLEL=$1
	VCF_FILE=$2
	Filter=$(echo "$3")
	VCF_OUT=$4
	DataName=$5
	CutoffCode=$6
	NumProc=$7
	VCFFIXUP=$8

	if [ $PARALLEL == "FALSE" ]; then 
		if [  $VCFFIXUP == "TRUE" ]; then vcffilter -s -f "$Filter" $VCF_FILE | vcffixup - > $VCF_OUT ; fi
		if [  $VCFFIXUP == "FALSE" ]; then vcffilter -s -f "$Filter" $VCF_FILE > $VCF_OUT; fi
		mawk '!/#/' $VCF_OUT | wc -l
	else
		tabix -H $VCF_FILE > ${VCF_OUT%.*}.header.vcf
		NumHeaderLines=$(tabix -H $VCF_FILE | wc -l)
		NumHeaderLines=$(($NumHeaderLines+2))
		if [  $VCFFIXUP == "TRUE" ]; then ls $DataName.$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcffilter -s -f $Filter | vcffixup - | tail -n +$NumHeaderLines" | cat ${VCF_OUT%.*}.header.vcf - | bgzip -@ $NumProc -c > ${VCF_OUT}.gz; fi
		if [  $VCFFIXUP == "FALSE" ]; then ls $DataName.$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcffilter -s -f $Filter | tail -n +$NumHeaderLines" | cat ${VCF_OUT%.*}.header.vcf - | bgzip -@ $NumProc -c > ${VCF_OUT}.gz; fi
		tabix -f -p vcf $VCF_OUT.gz
		#mawk '!/#/' $VCF_OUT.gz | wc -l
		rm ${VCF_OUT%.*}.header.vcf
	fi
}

###################################################################################################################
#stats functions
###################################################################################################################

function STATS() {
	name=$1[@]
	FILTERS=("${!name}")
	CutoffCode=$2
	BAM_PATH=$3
	VCF_FILE=$4
	REF_FILE=$5
	PopMap=$6
	CONFIG_FILE=$7
	HWE_SCRIPT=$8
	RADHAP_SCRIPT=$9
	DataName=${10}
	NumProc=${11}
	PARALLEL=${12}
	
	echo ${FILTERS[@]} | sed 's/ /\n/g' > $DataName.$CutoffCode.filters.txt
	ls $DataName.$CutoffCode.*vcf | parallel --no-notice -k -j $NumProc "grep -c '^dDocent_Contig' {} " | sort -r > $DataName.$CutoffCode.SNPcnt.txt
	ls $DataName.$CutoffCode.*vcf | parallel --no-notice -k -j $NumProc "grep '^dDocent_Contig' {} | cut -f1 | uniq | wc -l" | sort -r > $DataName.$CutoffCode.ContigCnt.txt

	#num missing genotypes for a particular snp
	zgrep -P 'dDocent_Contig_1\t' OpihiSK2014.A.25.10.Fltr18.HWE.recode.vcf.gz | sed '1q;d' | sed 's/\t/\n/g' | grep -c '^\.'
	
}


###################################################################################################################
#help info / manual
###################################################################################################################

NAME="$(basename "$0") v$VERSION  -- a program to filter vcf files with RAD data"

SYNOPSIS="$(basename "$0") [filter settings] [input files] [output file prefix] [parallelization]"

read -d '' DESCRIPTION <<"BLOCK"
fltrVCF is a tool to filter VCF files created by dDocentHPC. The filters can be run in any order.

        Arguments can be controlled from either the command line or a configuration file.  Filter
        thresholds can only be altered in the config.fltr file. Filters are described  and defined
		in the config.fltr file.

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
                -f <arg>        if set, controls filters to be run, in order. Argument should be 2
                                 digit numbers separated by spaces. -f "01 04 02"  or  -f 01\ 04\ 02
                                 will specify that filters 01, 04, and 02 will be run in succession.
                                 Filters are described in the config files. If -f is not set, the
                                 config file is used to determine the filters and order. If -f is
                                 set, it will override the config file. []
                -s <arg>        file with filter settings [config.fltr.clean.ind]

		[input files]
                -c <arg>        cutoff values used for reference genome [3.3]
                -b <arg>        path to mapping directory with *.bam and *.bed files [../mapping]
				-d <arg>		bed file describing complete data set. Required only if -P is set. 
								 [${b}/mapped.${c}.bed]
                -v <arg>        vcf file to be filtered [${b}/TotalRawSNPs.${c}.vcf]
                -g <arg>        reference genome fasta file [${b}/reference.${c}.fasta]
                -p <arg>        popmap file to use for defining population affiliation
                                 [${b}/popmap.${c}]
                -w <arg>        filter_hwe perl script [filter_hwe_by_pop_HPC.pl]
                -r <arg>        rad_haplotyper perl script [rad_haplotyperHPC116.pl]

        [output file prefix]
                -o <arg>        optional, all output files will be prefixed with this argument []

        [parallelization]
                -P				run every filter in parallel using GNU parallel. Requires *.bed files.  
								 If not set, then only natively-parallel filters will use multiple 
								 threads if -t > 1. Requires -d. [not set]
				-t <arg>        number of threads available for parallel processing [1]

EXAMPLES
		The following command is recommended for most users
				fltrVCF.bash -P -s config.fltr.ind
				
        The following two commands are the same, the first takes advantage of the defaults,
        the second does not.

                fltrVCF.bash -f "01 02 03" -c 25.10 -o ProjectX.A -t 40

                fltrVCF.bash -f "01 02 03" -c 25.10 -m ../mapping -v ../mapping/TotalRawSNPs.3.6.vcf
                        -p ../mapping/popmap.25.10 -s config.fltr.clean -w filter_hwe_by_pop.pl
                        -r rad_haplotyperHPC116.pl -o ProjectX.A -t 40

BLOCK


###################################################################################################################
#Read in command line arguments
###################################################################################################################
echo ""; echo $NAME; echo ""
echo "Dependencies:"
echo "	vcftools"
echo "	vcflib"
echo "	samtools"
echo "	perl"
echo "	mawk"
echo "	parallel"
echo "	rad_haplotyper116HPC" rad
echo "	filter_hwe_by_pop_HPC"
echo ""

echo "Reading options from command line:"
while getopts ":f:c:b:d:v:g:p:s:w:r:o:t:Ph" opt; do
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
    f)
        echo "	Filter order:                 $OPTARG"
        FILTERS=($OPTARG)
        ;;
	\?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    c)
        echo "" >&2
        echo "	Cutoffs:                  $OPTARG"
        CutoffCode=$OPTARG >&2
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    v)
        echo "	VCF File:                 $OPTARG"
        VCF=$OPTARG
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
    g)
        echo "	REF File:                 $OPTARG"
        REF_FILE=$OPTARG
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-(                 Option -$OPTARG requires an argument." >&2
        ;;
	b)
        echo "	Path to BAM files:        $OPTARG"
        BAM_PATH=$OPTARG
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    d)
        echo "	PopMap File:              $OPTARG"
        BED_FILE=$OPTARG
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    p)
        echo "	PopMap File:              $OPTARG"
        PopMap=$OPTARG
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    s)
        echo ""
		echo "	Settings File:            $OPTARG"
        CONFIG_FILE=$OPTARG
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    w)
        echo "	HWE Script:               $OPTARG"
        HWE_SCRIPT=$(${OPTARG})
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    r)
        echo "	Rad_Haplotyper script:    $OPTARG"
        RADHAP_SCRIPT=$(${OPTARG})
        ;;
    \?)
        echo "	ERROR :-/                 Invalid option: -$OPTARG" >&2
        ;;
    :)
        echo "	ERROR :-/                 Option -$OPTARG requires an argument." >&2
        ;;
    o)
        echo "	Output file prefix:       $OPTARG"
        DataName=$OPTARG
        ;;
    \?)
      echo "	ERROR :-/                   Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "	ERROR :-/                   Option -$OPTARG requires an argument." >&2
      ;;
    t)
      echo "	Number of threads:          $OPTARG"
      ;;
    \?)
      echo "	ERROR :-/                   Invalid option: -$OPTARG" >&2
      ;;
    :)
      echo "	ERROR :-/                   Option -$OPTARG requires an argument." >&2
      ;;
	P)
        echo "	Forcing all filters to run in parallel:	$OPTARG"
		PARALLEL=TRUE
		;;
  esac
done

echo "" >&2


###################################################################################################################
#Use the user input, config file, and defaults to control settings
###################################################################################################################

echo "Reading options from config file and setting defaults"
if [ -z ${CONFIG_FILE+x} ]; then 
	CONFIG_FILE=$(ls config.fltr.clean.ind)
	echo "	Settings are being loaded from: '${CONFIG_FILE}' by default"
else 
	CONFIG_FILE=$(ls ${CONFIG_FILE}) 
	echo "	SETTINGS are being loaded from file: '${CONFIG_FILE}'"
fi
if [ $CONFIG_FILE == "" ]; then 
	echo "	ERROR :-< 	configuration file is missing" >&2
	exit
fi

if [ -z ${FILTERS+x} ]; then 
	FILTERS=($(grep -P '^\t* *fltrVCF\t* *-f\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	echo "	Filters are set to '${FILTERS[@]}'"
else 
	echo "	FILTERS are set to '${FILTERS[@]}'"
fi
if [ ${FILTERS[0]} == "" ]; then 
	echo "	ERROR :-< 	filter settings (-f) are missing.  Please modify -f argument in config file or at the command line" >&2
	exit
fi

if [ -z ${CutoffCode+x} ]; then 
	CutoffCode=($(grep -P '^\t* *fltrVCF\t* *-c\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $CutoffCode == "" ]; then
		CutoffCode="3.3" 
	fi
fi
echo "	CutoffCode is set to '${CutoffCode}'"

if [ -z ${BAM_PATH+x} ]; then 
	BAM_PATH=($(grep -P '^\t* *fltrVCF\t* *-b\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $BAM_PATH == "" ]; then
		BAM_PATH="../mapping" 
	fi
fi
echo "	BAM_PATH is set to '${BAM_PATH}'"

if [ -z ${VCF_FILE+x} ]; then 
	VCF_FILE=($(grep -P '^\t* *fltrVCF\t* *-v\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $VCF_FILE == "" ]; then
		VCF_FILE=${BAM_PATH}"/TotalRawSNPs."${CutoffCode}".vcf"
	fi
	if [ ! -f $VCF_FILE ]; then
		echo "	ERROR:-<	$VCF_FILE does not exist"
		exit
	fi
fi
echo "	VCF_FILE is set to '${VCF_FILE}'"

if [ -z ${BED_FILE+x} ]; then 
	BED_FILE=($(grep -P '^\t* *fltrVCF\t* *-d\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ -z ${BED_FILE} ]; then
		BED_FILE=${BAM_PATH}"/mapped."${CutoffCode}".bed" 
	fi
fi	
echo "	Bed file is set to '${BED_FILE}'"

if [ -z ${REF_FILE+x} ]; then 
	REF_FILE=($(grep -P '^\t* *fltrVCF\t* *-g\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $REF_FILE == "" ]; then
		REF_FILE=${BAM_PATH}"/reference."${CutoffCode}".fasta" 
	fi
fi	
echo "	Reference genome is set to '${REF_FILE}'"

if [ -z ${PopMap+x} ]; then 
	PopMap=($(grep -P '^\t* *fltrVCF\t* *-p\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $PopMap == "" ]; then
		PopMap=${BAM_PATH}"/popmap."${CutoffCode} 
	fi
fi
echo "	PopMap is set to '${PopMap}'"

if [ -z ${HWE_SCRIPT+x} ]; then 
	HWE_SCRIPT=($(grep -P '^\t* *fltrVCF\t* *-w\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $HWE_SCRIPT == "" ]; then
		HWE_SCRIPT="filter_hwe_by_pop_HPC.pl" 
	fi
fi
echo "	HWE_SCRIPT is set to '${HWE_SCRIPT}'"

if [ -z ${RADHAP_SCRIPT+x} ]; then 
	RADHAP_SCRIPT=($(grep -P '^\t* *fltrVCF\t* *-r\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $RADHAP_SCRIPT == "" ]; then
		RADHAP_SCRIPT="rad_haplotyper116HPC.pl" 
	fi
fi
echo "	RADHAP_SCRIPT is set to '${RADHAP_SCRIPT}'"

if [ -z ${DataName+x} ]; then 
	DataName=($(grep -P '^\t* *fltrVCF\t* *-o\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $DataName == "" ]; then
		echo "	No output file prefix is specified.'" 
	else
		echo "	Output file prefix is set to '${DataName}'"
	fi
else
	echo "	Output file prefix is set to '${DataName}'"
fi


if [ -z ${NumProc+x} ]; then 
	NumProc=($(grep -P '^\tfltrVCF -t ' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [ $NumProc == "" ]; then
		NumProc=1 
	fi
fi
echo "	The number of threads is set to '${NumProc}'"

if [ ${FILTERS[0]} == "" ]; then 
	echo "	ERROR :-< 	filter settings (-f) are missing.  Please modify -f argument in config file or at the command line" >&2
	exit
fi

if [ ${PARALLEL} == "TRUE" ]; then
	if [ ! -f ${BED_FILE} ]; then
		echo "	ERROR :-\						Could not find *.bed file. Check the following settings -c -b -d -P."
		echo "									For help, try $(basename "$0") -h"
		exit
	fi
	#split the bed file by the number of processors, making sure that no contig is split between the files
	NumBedLines=$(($(wc -l ${BED_FILE} | sed 's/ .*//g')/${NumProc}))
	RemainderBedLines=$((${NumBedLines}%2))
	if [ ${RemainderBedLines} != 0 ]; then NumBedLines=$((NumBedLines+1)); fi
	split --lines=${NumBedLines} --numeric-suffixes --suffix-length=4 ${BED_FILE} $DataName.$CutoffCode.
	ls $DataName.$CutoffCode.[0-9][0-9][0-9][0-9] | parallel --no-notice -j $NumProc mv {} {}.bed
	
	echo "	-P set. All filters will be coerced to run in parallel."
	if [ "${VCF_FILE##*.}" == "vcf" ]; then
		if [ ! -f "${VCF_FILE}.gz" ]; then
			echo ""; echo "	" `date` " Using bgzip and tabix to compress and index VCF for parallelization"
			bgzip -@ $NumProc -c ${VCF_FILE} > ${VCF_FILE}.gz
			VCF_FILE=${VCF_FILE}.gz 
			tabix -p vcf ${VCF_FILE}
		else
			echo "	Changing -v to existing file ${VCF_FILE}.gz for parallelization.  If you didn't want this to happen, turn off -P or change -v"
			VCF_FILE=${VCF_FILE}.gz
		fi
	elif [ "${VCF_FILE##*.}" == "gz" ]; then
		if [ ! -f "${VCF_FILE}" ]; then
			echo "	ERROR:-(   this should be impossible to trigger"
		else
			echo "	${VCF_FILE} will be used for parallelization."
		fi
		if [ ! -f "${VCF_FILE}.tbi" ]; then
			echo "	${VCF_FILE}.tbi not found. Using tabix to index the *vcf.gz file"
			tabix -p vcf ${VCF_FILE}
		else
			echo "	${VCF_FILE}.tbi was successfully located."
		fi
	fi
else 
	echo "	-P not set by user. Only filters that natively support parallelization will be run in parallel."
	PARALLEL=FALSE
fi

echo ""

###################################################################################################################
#Run script
###################################################################################################################

MAIN FILTERS $CutoffCode $BAM_PATH $VCF_FILE $REF_FILE $PopMap $CONFIG_FILE $HWE_SCRIPT $RADHAP_SCRIPT $DataName $NumProc $PARALLEL

#cleanup files
ls $DataName.$CutoffCode.[0-9][0-9][0-9][0-9].bed | parallel --no-notice -j $NumProc "rm {}"
echo ""; echo `date` " --------------------------- Filtering complete! ---------------------------"; echo "" 