#!/bin/bash

#Improvements yet to be made
#	Shuffle the bed file to improve performance
#	Finish parallel option for rad_haplotyper
#	Enable option to make haplotype files with radhaplotyper
#	figure out how to combine haplotype files made by radhaplotyper in parallel mode

VERSION=4.2
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
	echo `date` " reading variables into MAIN"
	# name=$1[@]
	# FILTERS=("${!name}")
	# CutoffCode=$2
	# BAM_PATH=$3
	# VCF_FILE=$4
	# REF_FILE=$5
	# PopMap=$6
	# CONFIG_FILE=$7
	# HWE_SCRIPT=$8
	# RADHAP_SCRIPT=$9
	# DataName=${10}
	# NumProc=${11}
	# PARALLEL=${12}
	
	# echo "     Main filters are: ${FILTERS[@]}"
	# echo "          $CutoffCode"
	# echo "          $BAM_PATH"
	# echo "          $VCF_FILE"
	# echo "          $REF_FILE"
	# echo "          $PopMap"
	# echo "          $CONFIG_FILE"
	# echo "          $HWE_SCRIPT"
	# echo "          $RADHAP_SCRIPT"
	# echo "          $DataName"
	# echo "          $NumProc"
	# echo "          $PARALLEL"
	
	# if [[ $PARALLEL == "FALSE" ]]; then
		# VCF_FILE=$(ls -t ${DataName}*vcf | head -n 1)
		# VCF_FILE_2=$(ls -t ${DataName}*vcf | sed -n 2p)
	# else
		# VCF_FILE=$(ls -t ${DataName}*vcf.gz | head -n 1)
		# VCF_FILE_2=$(ls -t ${DataName}*vcf.gz | sed -n 2p)
	# fi
	
	#this keeps track of number of filters run
	COUNTER=1
	
	#this keeps track of number of times each filter is run
	COUNTERS=()
	for i in $(seq 1 100); do
		COUNTERS+=(0)
	done
	for i in ${FILTERS[@]}; do
		#echo $i Begin
		FILTER $i #$CutoffCode $BAM_PATH $VCF_FILE $REF_FILE $PopMap $CONFIG_FILE $HWE_SCRIPT $RADHAP_SCRIPT $DataName $NumProc $PARALLEL $VCF_FILE_2 
		#echo $i complete
		if [[ $PARALLEL == "FALSE" ]]; then
			VCF_FILE=$(ls -t ${DataName}*vcf | head -n 1)
			VCF_FILE_2=$(ls -t ${DataName}*vcf | sed -n 2p)
		else
			VCF_FILE=$(ls -t ${DataName}*vcf.gz | head -n 1)
			VCF_FILE_2=$(ls -t ${DataName}*vcf.gz | sed -n 2p)
		fi
		COUNTER=$((COUNTER+1))
	done

}


###################################################################################################################
#Define filters 
###################################################################################################################
function PARSE_THRESHOLDS(){
	#the purpose of this function is to parse multiple settings for 1 filter run several times
	local THRESH=$1
	if [[ $THRESH == *":"* ]]; then
		echo $(echo $THRESH | tr ':' '\t' | cut -f${COUNTERS[$INDEX]})
	else
		echo $THRESH
	fi
}


function FILTER(){
	FILTER_ID=$1
	# CutoffCode=$2
	# BAM_PATH=$3
	# VCF_FILE=$4
	# REF_FILE=$5
	# PopMap=$6
	# CONFIG_FILE=$7
	# HWE_SCRIPT=$8
	# RADHAP_SCRIPT=$9
	# DataName=${10}
	# NumProc=${11}
	# PARALLEL=${12}
	# VCF_FILE_2=${13}
	
	# echo "     $FILTER_ID"
	# echo "     $CutoffCode"
	# echo "     $BAM_PATH"
	# echo "     $VCF_FILE"
	# echo "     $REF_FILE"
	# echo "     $PopMap"
	# echo "     $CONFIG_FILE"
	# echo "     $HWE_SCRIPT"
	# echo "     $RADHAP_SCRIPT"
	# echo "     $DataName"
	# echo "     $NumProc"
	# echo "     $PARALLEL"
	# echo "     $COUNTER"
	
	#Keep track of times each filter has been run
	INDEX=$(echo $FILTER_ID | sed -e 's:^0*::')
	COUNTERS[$INDEX]=$((1+COUNTERS[$INDEX]))

	
	VCF_OUT=$DataName$CutoffCode.Fltr${FILTER_ID}.${COUNTER}
	if [[ $FILTER_ID == "01" ]]; then
		echo; echo `date` "---------------------------FILTER01: Remove sites with Y < alleles < X -----------------------------"
		COUNTER01=$((1+COUNTER01))
		#get settings from config file
		THRESHOLD=$(grep -P '^\t* *01\t* *vcftools\t* *--min-alleles' ${CONFIG_FILE} | sed 's/\t* *01\t* *vcftools\t* *--min-alleles\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLD" ]]; then THRESHOLD=2; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		THRESHOLDb=$(grep -P '^\t* *01\t* *vcftools\t* *--max-alleles' ${CONFIG_FILE} | sed 's/\t* *01\t* *vcftools\t* *--max-alleles\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLDb" ]]; then THRESHOLDb=2; fi
		THRESHOLDb=$(PARSE_THRESHOLDS $THRESHOLDb) 
		#Filter="\"--min-alleles $THRESHOLD --max-alleles $THRESHOLDb --recode --recode-INFO-all\""
		Filter="--min-alleles $THRESHOLD --max-alleles $THRESHOLDb --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr${FILTER_ID}.${COUNTER}
		#call function to filter vcf
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc

	elif [[ $FILTER_ID == "02" ]]; then
		echo; echo `date` "---------------------------FILTER02: Remove Sites with Indels -----------------------------"
		Filter="--remove-indels --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "03" ]]; then
		echo; echo `date` "---------------------------FILTER03: Remove Sites with QUAL < minQ -----------------------------"
		THRESHOLD=($(grep -P '^\t* *03\t* *vcftools\t* *--minQ' ${CONFIG_FILE} | sed 's/\t* *03\t* *vcftools\t* *--minQ\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=30; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="--minQ ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "04" ]]; then
		echo; echo `date` "---------------------------FILTER04: Remove Sites With Mean Depth of Coverage < min-meanDP -----------------------------"
		THRESHOLD=($(grep -P '^\t* *04\t* *vcftools\t* *--min-meanDP' ${CONFIG_FILE} | sed 's/\t* *04\t* *vcftools\t* *--min-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=2; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="--min-meanDP ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "05" ]]; then
		echo; echo `date` "---------------------------FILTER05: Remove sites called in <X proportion of individuals -----------------------------"
		THRESHOLD=($(grep -P '^\t* *05\t* *vcftools\t* *--max-missing' ${CONFIG_FILE} | sed 's/\t* *05\t* *vcftools\t* *--max-missing\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.5; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="--max-missing ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "06" ]]; then
		echo; echo `date` "---------------------------FILTER06: Remove sites with Average Allele Balance deviating too far from 0.5 while keeping those with AB=0  -----------------------------"
		THRESHOLD=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *min' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.375; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		THRESHOLDb=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *max' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *max\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDb}" ]]; then ${THRESHOLDb}=0.625; fi
		THRESHOLDb=$(PARSE_THRESHOLDS $THRESHOLDb) 
		Filter="AB > $THRESHOLD & AB < $THRESHOLDb | AB = 0"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER "TRUE" #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "TRUE" #the last option "TRUE" is for vcffixup

#UNDER CONSTRUCTION V
#	elif [[ $FILTER_ID == "96" ]]; then
#		echo; echo `date` "---------------------------FILTER06: Remove contigs with Average Allele Balance deviating too far from 0.5 while keeping those with AB=0  -----------------------------"
#		THRESHOLD=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *min' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
#		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.375; fi
#		THRESHOLDb=($(grep -P '^\t* *06\t* *vcffilter\t* *AB\t* *max' ${CONFIG_FILE} | sed 's/\t* *06\t* *vcffilter\t* *AB\t* *max\t* *//g' | sed 's/\t* *#.*//g' )) 
#		if [[ -z "${THRESHOLDb}" ]]; then ${THRESHOLDb}=0.625; fi
#		Filter="\"AB > $THRESHOLD & AB < $THRESHOLDb | AB = 0\""
#		VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
#		#get the AB data
#		grep '^dDocent_Contig' $VCF_FILE | cut -f1,8 | sed 's/;/\t/g' | cut -f1,2 | sed 's/AB=//' > ${VCF_OUT%.*}.AB.txt
#		#calculate mean AB
#		awk '{print $1}' AB.txt | uniq | uniq.contigs.txt
#		awk '$1 == "dDocent_Contig_1"' AB.txt | awk '$2 != 0' | awk '{total += $2; n++ } END {if (n > 0) print "dDocent_Contig_1 " total / n} '
#		#do it in parallel
#		#cat uniq.contigs.txt | parallel --no-notice -k -j 
#		FILTER_VCFFILTER $PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "TRUE" #the last option "TRUE" is for vcffixup
#UNDER CONSTRUCTION ^

	elif [[ $FILTER_ID == "07" ]]; then
		echo; echo `date` "---------------------------FILTER07: Remove sites with Alternate Allele Count <=X -----------------------------"
		THRESHOLD=($(grep -P '^\t* *07\t* *vcffilter\t* *AC\t* *min' ${CONFIG_FILE} | sed 's/\t* *07\t* *vcffilter\t* *AC\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=1; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="AC > $THRESHOLD & AN - AC > $THRESHOLD"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER "FALSE" #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [[ $FILTER_ID == "08" ]]; then
		echo; echo `date` "---------------------------FILTER08: Remove sites covered by both F and R reads-----------------------------"
		# This should not be run if you expect to have a lot of overlap between forward and rev reads (Miseq 2 x 300bp)
		THRESHOLD=($(grep -P '^\t* *08\t* *vcffilter\t* *SAF\/SAR\t* *min' ${CONFIG_FILE} | sed 's/\t* *08\t* *vcffilter\t* *SAF\/SAR\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=1; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="SAF / SAR > ${THRESHOLD} & SRF / SRR > ${THRESHOLD} | SAR / SAF > ${THRESHOLD} & SRR / SRF > ${THRESHOLD}"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER "FALSE" #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [[ $FILTER_ID == "09" ]]; then
		echo; echo `date` "---------------------------FILTER09: Remove sites with low/high ratio of mean mapping quality of alt to ref -----------------------------"
		# The next filter looks at the ratio of mapping qualities between reference and alternate alleles.  The rationale here is that, again, because RADseq 
		# loci and alleles all should start from the same genomic location there should not be large discrepancy between the mapping qualities of two alleles.
		THRESHOLD=($(grep -P '^\t* *09\t* *vcffilter\t* *MQM\/MQMR\t* *min' ${CONFIG_FILE} | sed 's/\t* *09\t* *vcffilter\t* *MQM\/MQMR\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.1; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		calc(){ awk "BEGIN { print "$*" }"; }
		THRESHOLDmin=$(calc 1-${THRESHOLD})
		#THRESHOLDmin=$(echo 1-${THRESHOLD} | R --vanilla --quiet | sed -n '2s/.* //p')
		THRESHOLDmax=$(calc 1/${THRESHOLDmin})
		#THRESHOLDmax=$(echo 1/${THRESHOLDmin} | R --vanilla --quiet | sed -n '2s/.* //p')
#		Filter="\"MQM / MQMR > 1 - ${THRESHOLDmin} & MQM / MQMR < ${THRESHOLDmax}\""
		Filter="MQM / MQMR > ${THRESHOLDmin} & MQM / MQMR < ${THRESHOLDmax}"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER "FALSE" #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [[ $FILTER_ID == "10" ]]; then
		echo; echo `date` "---------------------------FILTER10: Remove sites with one allele only supported by reads not properly paired -----------------------------"
		# another filter that can be applied is whether or not their is a discrepancy in the properly paired status  for reads supporting reference or alternate alleles.
		# Since de novo assembly is not perfect, some loci will only have unpaired reads mapping to them. This is not a problem. The problem occurs when all the reads supporting 
		# the reference allele are paired but not supporting the alternate allele. That is indicative of a problem.
		Filter="PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER "FALSE" #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"
		
	elif [[ $FILTER_ID == "11" ]]; then
		echo; echo `date` "---------------------------FILTER11: Remove sites with Quality/DP ratio < 0.25 -----------------------------"
		# There is a positive relationship between depth of coverage and inflation of quality score. JPuritz controls for this in two ways.
		# first, by removing any locus that has a quality score below 1/4 of the read depth.
		THRESHOLD=($(grep -P '^\t* *11\t* *vcffilter\t* *QUAL\/DP\t* *min' ${CONFIG_FILE} | sed 's/\t* *11\t* *vcffilter\t* *QUAL\/DP\t* *min\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.25; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="QUAL / DP > ${THRESHOLD}"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID.vcf
		FILTER_VCFFILTER "FALSE" #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc "FALSE"

	elif [[ $FILTER_ID == "12" ]]; then
		echo; echo `date` "---------------------------FILTER12: Remove sites with Quality > 2xDP -----------------------------"
		# second, is a multistep process
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		Filter="--exclude-positions $VCF_OUT.lowQDloci --recode --recode-INFO-all"
		# create a list of the depth of each locus
		if [[ $PARALLEL == "FALSE" ]]; then
			cut -f8 $VCF_FILE | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $VCF_OUT.DEPTH
		else
			zcat $VCF_FILE | grep -oe "DP=[0-9]*" | sed -s 's/DP=//g' > $VCF_OUT.DEPTH
		fi
		#Then make a list of quality scores
		mawk '!/#/' $VCF_FILE | cut -f1,2,6 > $VCF_OUT.loci.qual
		#Then calculate mean depth
		MeanDepth=$(mawk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $VCF_OUT.DEPTH )
		#echo $MeanDepth
		#MeanDepth2=$(python -c "print int($MeanDepth+3*($MeanDepth**0.5))" )
		MeanDepth2=$(echo "$MeanDepth+3*sqrt($MeanDepth)" | bc)
		#echo $MeanDepth2
		#Next we paste the depth and quality files together and find the loci above the cutoff that do not have quality scores 2 times the depth
		paste $VCF_OUT.loci.qual $VCF_OUT.DEPTH | mawk -v x=$MeanDepth2 '$4 > x' | mawk '$3 < 2 * $4' > $VCF_OUT.lowQDloci
		#cat $VCF_OUT.lowQDloci
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "13" ]]; then
		echo; echo `date` "---------------------------FILTER13: Remove sites with mean DP greater than X  -----------------------------"
		THRESHOLD=($(grep -P '^\t* *13\t* *vcftools\t* *--max-meanDP' ${CONFIG_FILE} | sed 's/\t* *13\t* *vcftools\t* *--max-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=250; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter2="--site-depth"
		Filter="--max-meanDP $THRESHOLD --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		#calculate the depth across loci with VCFtools
		if [[ $PARALLEL == "FALSE" ]]; then 
			vcftools --vcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			#Now let’s take VCFtools output and cut it to only the depth scores
			cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			#Now let’s calculate the average depth by dividing the above file by the number of individuals
			NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $VCF_FILE )
		else
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			#Now let’s take VCFtools output and cut it to only the depth scores
			cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			#Now let’s calculate the average depth by dividing the above file by the number of individuals
			NumInd=$(gunzip -c $VCF_FILE | awk '{if ($1 == "#CHROM"){print NF-9; exit}}' )
		fi
		mawk '!/D/' $VCF_OUT.site.depth | mawk -v x=$NumInd '{print $1/x}' > $VCF_OUT.meandepthpersitebefore
		cp $VCF_OUT.meandepthpersitebefore meandepthpersitebefore
		
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [0:*] 
unset label
set title "Histogram of mean depth per site before filter"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 5
plot 'meandepthpersitebefore' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

   
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

		if [[ $PARALLEL == "FALSE" ]]; then 
			vcftools --vcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			#Now let’s take VCFtools output and cut it to only the depth scores
			cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			#Now let’s calculate the average depth by dividing the above file by the number of individuals
			NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $VCF_FILE )
		else
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			#Now let’s take VCFtools output and cut it to only the depth scores
			cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			#Now let’s calculate the average depth by dividing the above file by the number of individuals
			NumInd=$(gunzip -c $VCF_FILE | awk '{if ($1 == "#CHROM"){print NF-9; exit}}' )
		fi
		mawk '!/D/' $VCF_OUT.site.depth | mawk -v x=$NumInd '{print $1/x}' > $VCF_OUT.meandepthpersiteafter
		cp $VCF_OUT.meandepthpersiteafter meandepthpersiteafter
		
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [0:*] 
unset label
set title "Histogram of mean depth per site after filter"
set ylabel "Number of Occurrences"
set xlabel "Mean Depth"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 5
plot 'meandepthpersiteafter' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
		
		rm meandepthpersite*
		
		
	elif [[ $FILTER_ID == "14" ]]; then
		echo; echo `date` "---------------------------FILTER14: If individual's genotype has DP < X, convert to missing data -----------------------------"
		THRESHOLD=($(grep -P '^\t* *14\t* *vcftools\t* *--minDP' ${CONFIG_FILE} | sed 's/\t* *14\t* *vcftools\t* *--minDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=3; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="--minDP ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "15" ]]; then
		echo; echo `date` "---------------------------FILTER15: Remove sites with maf > minor allele frequency > max-maf -----------------------------"
		THRESHOLDa=($(grep -P '^\t* *15\t* *vcftools\t* *--maf' ${CONFIG_FILE} | sed 's/\t* *15\t* *vcftools\t* *--maf\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDa}" ]]; then ${THRESHOLDa}=0.005; fi
		THRESHOLDa=$(PARSE_THRESHOLDS $THRESHOLDa) 
		THRESHOLDb=($(grep -P '^\t* *15\t* *vcftools\t* *--max-maf' ${CONFIG_FILE} | sed 's/\t* *15\t* *vcftools\t* *--max-maf\t* *//g' | sed 's/\t* *#.*//g' )) 
		THRESHOLDb=$(PARSE_THRESHOLDS $THRESHOLDb) 
		if [[ -z "${THRESHOLDb}" ]]; then ${THRESHOLDb}=0.995; fi
		# Remove sites with minor allele frequency: maf < x < max-maf
		# inspect the AF values in the vcf.  This will affect the frequency of rare variants
		Filter="--maf ${THRESHOLDa} --max-maf ${THRESHOLDb} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "16" ]]; then
		echo; echo `date` "---------------------------FILTER16: Remove Individuals with Too Much Missing Data-----------------------------"
		THRESHOLD=($(grep -P '^\t* *16\t* *vcftools\t* *--missing-indv' ${CONFIG_FILE} | sed 's/\t* *16\t* *vcftools\t* *--missing-indv\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.5; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter="--missing-indv"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		if [[ $PARALLEL == "FALSE" ]]; then 
			vcftools --vcf ${VCF_FILE} $Filter --out $VCF_OUT 2> /dev/null
		else
			vcftools --gzvcf ${VCF_FILE} $Filter --out $VCF_OUT 2> /dev/null
		fi
		#echo; echo "  Missing Data Report, file=*out.imiss, Numbers Near 0 are Good"
		#cat $VCF_OUT.imiss
		#graph missing data for individuals
		mawk '!/IN/' $VCF_OUT.imiss | cut -f5 | sort -r > $VCF_OUT.totalmissing
		cp $VCF_OUT.totalmissing totalmissing
		seq 1 $(cat totalmissing | wc -l) > ind.txt
		paste ind.txt totalmissing > imiss.dat
gnuplot << \EOF 
		#filename=system("echo $VCF_OUT.totalmissing")
		set terminal dumb size 120, 30
		set autoscale 
		unset label
		set title "Histogram of % missing data per individual before filter. Bars to the left are desireable."
		set ylabel "Number of Individuals"
		set xlabel "% missing genotypes"
		set yrange [0:*]
		set xrange [0:1]
		binwidth=0.01
		bin(x,width)=width*floor(x/width) + binwidth/2.0
		plot 'totalmissing' using (bin($1,binwidth)):(1.0) smooth freq with boxes
		pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of % missing data per individual."
set ylabel "% of missing data"
set xlabel "Individual"
xmax="`cut -f1 imiss.dat | tail -1`"
xmax=xmax+1
set xrange [0:xmax]
set yrange [0:1]
plot 'imiss.dat' pt "*" 
pause -1
EOF
		##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#id individuals to remove based upon the proportion of missing data
		#examine the plot and the $DataName$CutoffCode.out.imiss file to determine if this cutoff is ok
		mawk -v x="$THRESHOLD" ' $5 > x' $VCF_OUT.imiss | cut -f1 > $VCF_OUT.lowDP-2.indv
		##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		Filter="--remove $VCF_OUT.lowDP-2.indv --recode --recode-INFO-all"
		#list of individuals to remove
		echo " Individuals with too much missing data:"
		cat $VCF_OUT.lowDP-2.indv
		#remove individuals with low reads
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 
		#need to remove individuals from header
		if [[ $PARALLEL == "TRUE" ]]; then gunzip $VCF_OUT.recode.vcf.gz; fi
		while read i; do
			sed -i "s/\t$i//g" $VCF_OUT.recode.vcf
		done < $VCF_OUT.lowDP-2.indv
		if [[ $PARALLEL == "TRUE" ]]; then 
			bgzip -@ $NumProc -c $VCF_OUT.recode.vcf > $VCF_OUT.recode.vcf.gz
			tabix -p vcf $VCF_OUT.recode.vcf.gz
			rm $VCF_OUT.recode.vcf
		fi
		rm totalmissing ind.txt imiss.dat
		
		
	elif [[ $FILTER_ID == "17" ]]; then
		echo; echo `date` "---------------------------FILTER17: Remove sites with data missing for too many individuals in a population -----------------------------"
		THRESHOLD=($(grep -P '^\t* *17\t* *vcftools\t* *--missing-sites' ${CONFIG_FILE} | sed 's/\t* *17\t* *vcftools\t* *--missing-sites\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.5; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		# restrict the data to loci with a high percentage of individuals that geneotyped
		#SitesMissData: on/off switch for this filter at beginning of script
		echo "     Using PopMap File: $PopMap"
		#cat $PopMap
		# get the popnames from the popmap
		popnames=`mawk '{print $2}' $PopMap | sort | uniq `
		missingDataByPop() {
			PARALLEL=$3
			VCF_FILE=$4
			VCF_OUT=$5
			mawk -v pop=$1 '$2 == pop' $2 > $VCF_OUT.$1.keep
			if [[ $PARALLEL == "FALSE" ]]; then
				vcftools --vcf $VCF_FILE --keep $VCF_OUT.$1.keep --missing-site --out $VCF_OUT.$1 2> /dev/null
			else
				vcftools --gzvcf $VCF_FILE --keep $VCF_OUT.$1.keep --missing-site --out $VCF_OUT.$1 2> /dev/null
			fi
			
#			cat $VCF_OUT.$1.lmiss | cut -f6 | tail -n +2 | sort -r > lmiss.txt
		}
		export -f missingDataByPop
		parallel --no-notice -k -j $NumProc "missingDataByPop {} $PopMap $PARALLEL $VCF_FILE $VCF_OUT 2> /dev/null" ::: ${popnames[*]}
		
		while re
		# cat opihiSK2014A.10.25.Fltr17.8.opihiSK2014-EastMaui1-A.lmiss | cut -f6 | tail -n +2 | sort -r > lmiss.txt
		# seq 1 $(cat lmiss.txt | wc -l) > snplmiss.txt
		# paste snplmiss.txt lmiss.txt > lmiss.dat
		
# gnuplot << \EOF 
# set terminal dumb size 120, 30
# set autoscale 
# unset label
# set title "Scatter plot of the proportion of missing genotypes by SNP before filtering."
# set ylabel "% missing genotypes"
# set xlabel "SNP"
# xmax="`cut -f1 lmiss.dat | tail -1`"
# xmax=xmax+1
# set xrange [0:xmax]
# set yrange [0:1]
# plot 'lmiss.dat' pt "*" 
# pause -1
# EOF
		
		# remove loci with more than X proportion of missing data
		# if you have mixed sequence lengths, this will affect if the longer regions are typed
		# UPDATE % missing threshold here	
		cat $VCF_OUT.*.lmiss | mawk '!/CHR/' | mawk -v x=$THRESHOLD '$6 > x' | cut -f1,2 | sort | uniq > $VCF_OUT.badloci
		
		Filter="--exclude-positions $VCF_OUT.badloci --recode --recode-INFO-all"
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 
		
		cut -f1 $VCF_OUT.badloci | parallel "grep -v {} "

	elif [[ $FILTER_ID == "18" ]]; then
		echo; echo `date` "---------------------------FILTER18: Remove sites not in HWE p<X) -----------------------------"
		THRESHOLD=($(grep -P '^\t* *18\t* *filter_hwe_by_pop_HPC' ${CONFIG_FILE} | sed 's/\t* *18\t* *filter_hwe_by_pop_HPC\t* *//g' | sed 's/\t* *#.*//g' )) 
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.001; fi
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		if [[ $PARALLEL == "TRUE" ]]; then 
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
		
		
	elif [[ $FILTER_ID == "19" ]]; then
		echo; echo `date` "---------------------------FILTER19: Run rad_haplotyper to id paralogs,create haplotypes, etc -----------------------------"
		THRESHOLDa=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-d\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-d\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDa}" ]]; then ${THRESHOLDa}=50; fi
		THRESHOLDa=$(PARSE_THRESHOLDS $THRESHOLDa) 
		THRESHOLDb=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-mp\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-mp\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDb}" ]]; then ${THRESHOLDb}=10; fi
		THRESHOLDb=$(PARSE_THRESHOLDS $THRESHOLDb) 
		THRESHOLDc=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-u\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-u\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDc}" ]]; then ${THRESHOLDd}=30; fi
		THRESHOLDc=$(PARSE_THRESHOLDS $THRESHOLDc) 
		THRESHOLDd=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-ml\t* *' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-ml\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDd}" ]]; then ${THRESHOLDd}=10; fi
		THRESHOLDd=$(PARSE_THRESHOLDS $THRESHOLDd) 
		THRESHOLDe=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-h\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-h\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDe}" ]]; then ${THRESHOLDe}=100; fi
		THRESHOLDe=$(PARSE_THRESHOLDS $THRESHOLDe) 
		THRESHOLDf=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-z\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-z\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDf}" ]]; then ${THRESHOLDf}=0.1; fi
		THRESHOLDf=$(PARSE_THRESHOLDS $THRESHOLDf) 
		THRESHOLDg=($(grep -P '^\t* *19\t* *rad_haplotyper\t* *-m\t* *\d' ${CONFIG_FILE} | sed 's/\t* *19\t* *rad_haplotyper\t* *-m\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLDg}" ]]; then ${THRESHOLDg}=0.5; fi
		THRESHOLDg=$(PARSE_THRESHOLDS $THRESHOLDg) 
		VCF_OUT=$DataName$CutoffCode
		echo "     Read Sampling Depth:					$THRESHOLDa"
		echo "     Max Paralogous Individuals: 				$THRESHOLDb"
		echo "     Max Num SNPs per Contig:				$THRESHOLDc"
		echo "     Max Low Cov Individuals:				$THRESHOLDd"
		echo "     Max Haplotypes in XS of NumSNPs:			$THRESHOLDe"
		echo "     Num or Prop of Reads With Rare Haplotypes to Ignore:	$THRESHOLDf"
		echo "     Min Prop of Individuals with Haplotypes to Keep Locus:	$THRESHOLDg"
		
		# restrict the data to loci with a high percentage of individuals that geneotyped
		###rad_haplotyper will create a vcf file with only the haplotypable snps, a genepop file with haps, and an ima file with haps
		#note, rad haplotyper skips complex polymorphisms by default
		#runs ceb version of rad haplotyper
		# if [[ $PARALLEL == "TRUE" ]]; then 
			# gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			# perl $RADHAP_SCRIPT -v ${VCF_FILE%.*} -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -co ${CutoffCode} -dn ${DataName}  -p ${PopMap} -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf -g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			# mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			# bgzip -@ $NumProc -c $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf > $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
			# tabix -p vcf $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
		# else
			# perl $RADHAP_SCRIPT -v $VCF_FILE -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -co ${CutoffCode} -dn ${DataName}  -p ${PopMap} -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf -g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			# mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
		# fi

		#runs official rad haplotyper with minor addition of a bam path -bp
		sed "s/\t/$CutoffCode\t/g" $PopMap > $PopMap$CutoffCode
		
		if [[ $PARALLEL == "TRUE" ]]; then 
			############## Works, but runs out of memory#######################
				# gunzip -c $VCF_FILE > ${VCF_FILE%.*}
				# #modify individual names in vcf to include cutoffcode
				# LineNum=$(awk '/^#CHROM/{print NR; exit}' ${VCF_FILE%.*})
				# FirstCols=$(grep '^#CHROM' ${VCF_FILE%.*} | cut -f1-9)
				# NewIndNames=$(grep '^#CHROM' ${VCF_FILE%.*} | cut -f10- | sed "s/\t/$CutoffCode\t/g" | sed "s/$/$CutoffCode/g")
				# NewLine=$(echo $FirstCols $NewIndNames)
				# sed "${LineNum}s/.*/$(echo $NewLine | sed "s/ /\t/g")/" ${VCF_FILE%.*} > ${VCF_FILE%.*}.1.vcf
				
				# perl $RADHAP_SCRIPT -v ${VCF_FILE%.*}.1.vcf -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}$CutoffCode -o ${VCF_OUT}.Fltr${FILTER_ID}.Haplotyped.vcf #-g ${VCF_OUT}.Fltr${FILTER_ID}.${PopMap##*/}.haps.genepop -a ${VCF_OUT}.Fltr${FILTER_ID}.${PopMap##*/}.haps.ima
				
				# mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
				# bgzip -@ $NumProc -c $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf > $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
				# tabix -p vcf $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz
				# rm ${VCF_FILE%.*}.1.vcf
			################ End ##################################################
			
			############# Alternative that splits up vcf using bed files, parallelized radhaplotyper, and compartmentalizes rad haplotyper file output in sub directories#########
				
				ForceParallelRadHap() {
					RADHAP_SCRIPT=$1
					VCF_FILE=$2
					NumProc=$3
					THRESHOLDa=$4
					THRESHOLDb=$5
					THRESHOLDc=$6
					THRESHOLDd=$7
					THRESHOLDe=$8
					THRESHOLDf=$9
					THRESHOLDg=${10}
				    REF_FILE=${11}
					BAM_PATH=${12}
					PopMap=${13}
					CutoffCode=${14}
					VCF_OUT=${15}
					FILTER_ID=${16}
					DataName=${17}
					Indexing=${18}
					
					echo $RADHAP_SCRIPT
					echo $VCF_FILE
					echo $NumProc
					echo $THRESHOLDa
					echo $THRESHOLDb
					echo $THRESHOLDc
					echo $THRESHOLDd
					echo $THRESHOLDe
					echo $THRESHOLDf
					echo $THRESHOLDg
					echo $REF_FILE
					echo $BAM_PATH
					echo $PopMap
					echo $CutoffCode
					echo $VCF_OUT
					echo $FILTER_ID
					echo $DataName
					echo $Indexing
 
					
					mkdir radhap${Indexing}
					cp $DataName$CutoffCode.${Indexing}.bed radhap${Indexing}
					cp $RADHAP_SCRIPT radhap${Indexing}
					cd radhap${Indexing}

					# tabix -H ../$VCF_FILE > $VCF_OUT.header.vcf
					# NumHeaderLines=$(tabix -H ../$VCF_FILE | wc -l)
					# NumHeaderLines=$(($NumHeaderLines+1))
					tabix -h -R ${DataName}${CutoffCode}.${Indexing}.bed ../$VCF_FILE > ${Indexing}.vcf
					
					#modify individual names in vcf to include cutoffcode
					LineNum=$(awk '/^#CHROM/{print NR; exit}' ${Indexing}.vcf)
					FirstCols=$(grep '^#CHROM' ${Indexing}.vcf | cut -f1-9)
					NewIndNames=$(grep '^#CHROM' ${Indexing}.vcf | cut -f10- | sed "s/\t/$CutoffCode\t/g" | sed "s/$/$CutoffCode/g")
					NewLine=$(echo $FirstCols $NewIndNames)
					sed -i "${LineNum}s/.*/$(echo $NewLine | sed "s/ /\t/g")/" ${Indexing}.vcf 
					
					#-x is the number of threads
					perl $RADHAP_SCRIPT -v ${Indexing}.vcf -x 1 -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}$CutoffCode -o ${VCF_OUT}.Fltr${FILTER_ID}.Haplotyped.vcf #-g ${VCF_OUT}.Fltr${FILTER_ID}.${PopMap##*/}.haps.genepop -a ${VCF_OUT}.Fltr${FILTER_ID}.${PopMap##*/}.haps.ima

					if [[ ! -f ../$VCF_OUT.Fltr$FILTER_ID.stats.out ]]; then head -n 2 stats.out > ../$VCF_OUT.Fltr$FILTER_ID.stats.out; fi
					tail -n +3 stats.out >> ../$VCF_OUT.Fltr$FILTER_ID.stats.out
					
					cd ..

				}
				export -f ForceParallelRadHap
				
				#make files to receive output from radhap instances
				touch $VCF_OUT.Fltr$FILTER_ID.stats.out
				
				seq -f "%04g" 0 $((NumProc-1)) | parallel --no-notice -k -j $NumProc "ForceParallelRadHap $RADHAP_SCRIPT $VCF_FILE $NumProc $THRESHOLDa $THRESHOLDb $THRESHOLDc $THRESHOLDd $THRESHOLDe $THRESHOLDf $THRESHOLDg $REF_FILE $BAM_PATH $PopMap $CutoffCode $VCF_OUT $FILTER_ID $DataName {} "
				
			#################END#######################################################
		else
			#modify individual names in vcf to include cutoffcode
			LineNum=$(awk '/^#CHROM/{print NR; exit}' $VCF_FILE)
			FirstCols=$(grep '^#CHROM' $VCF_FILE | cut -f1-9)
			NewIndNames=$(grep '^#CHROM' $VCF_FILE | cut -f10- | sed "s/\t/$CutoffCode\t/g" | sed "s/$/$CutoffCode/g")
			NewLine=$(echo $FirstCols $NewIndNames)
			sed "${LineNum}s/.*/$(echo $NewLine | sed "s/ /\t/g")/" $VCF_FILE > $VCF_FILE.1.vcf
			
			perl $RADHAP_SCRIPT -v $VCF_FILE.1.vcf -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}$CutoffCode -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf #-g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
			mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			rm $VCF_FILE.1.vcf
		fi
		
		rm $PopMap$CutoffCode

		#rename outputfiles
		if [[ -f "fail.log" ]]; then mv fail.log $VCF_OUT.Fltr$FILTER_ID.fail.log; fi
		if [[ -f "haplo_dump.out" ]]; then mv haplo_dump.out $VCF_OUT.Fltr$FILTER_ID.haplo_dump.out; fi
		if [[ -f "hap_log.out" ]]; then mv hap_log.out $VCF_OUT.Fltr$FILTER_ID.hap_log.out; fi
		if [[ -f "stats.out" ]]; then mv stats.out $VCF_OUT.Fltr$FILTER_ID.stats.out; fi
		if [[ -f "ind_stats.out" ]]; then mv ind_stats.out $VCF_OUT.Fltr$FILTER_ID.ind_stats.out; fi
		if [[ -f "snp_dump.out" ]]; then mv snp_dump.out $VCF_OUT.Fltr$FILTER_ID.snp_dump.out; fi
		if [[ -f "allele_dump.out" ]]; then mv allele_dump.out $VCF_OUT.Fltr$FILTER_ID.allele_dump.out; fi
		if [[ -f "temp.vcf" ]]; then mv temp.vcf $VCF_OUT.Fltr$FILTER_ID.temp.vcf; fi
		
		#Why and how many loci were removed by rad_haplotyper
		Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
		parallel --gnu --null "grep {} $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f1 > $VCF_OUT.Fltr$FILTER_ID.{}" ::: "${Explanations[@]}"
		parallel --gnu --null "echo -n {}' removed, ' && wc -l $VCF_OUT.Fltr$FILTER_ID.{}" ::: "${Explanations[@]}"
		echo file $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf has 
		mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
		echo SNPs

	elif [[ $FILTER_ID == "20" ]]; then
		echo; echo `date` "---------------------------FILTER20: Select 1 Random SNP per Contig -----------------------------"
		#Script to take random SNP from every contig in a vcffile
		#Get name of VCF file
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		if [[ $PARALLEL == "TRUE" ]]; then 
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
		if [[ $PARALLEL == "TRUE" ]]; then 
			bgzip -@ $NumProc -c $VCF_OUT.randSNPperLoc.vcf > $VCF_OUT.randSNPperLoc.vcf.gz
			tabix -p vcf $VCF_OUT.randSNPperLoc.vcf.gz
		fi

	elif [[ $FILTER_ID == "21" ]]; then
		echo; echo `date` "---------------------------FILTER21: Select Most Informative SNP per Contig -----------------------------"
		#Script to take random SNP from every contig in a vcffile
		#Get name of VCF file
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		if [[ $PARALLEL == "TRUE" ]]; then 
			gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			VCF_FILE=${VCF_FILE%.*}
		fi

		NAME=$(echo $VCF_FILE | sed -e 's/\.recode.*//g' | sed -e 's/.vcf//g' ) 
		#Calculate number of SNPs
		Loci=(`mawk '!/#/' $VCF_FILE | wc -l `)
		#Generate list of random numbers
		#seq 1 500000 | shuf | head -$Loci > $VCF_OUT.nq
		
		#Adjust AF values so that 0.5 is the max, 0.5 has most info
		AF=($(vcf-query *vcf -f '%INFO/AF\n'))
		AF_len=$(( ${#AF[@]} - 1))
		for i in $(seq 0 $AF_len); do
			if (( $(echo "${AF[$i]} > 0.5" | bc) )); then
				AF[$i]=$(echo "1 - ${AF[$i]}" | bc) 
			fi
		done
		echo ${AF[@]} | tr " " "\n" > $VCF_OUT.mi 
		
		
		#create temporary file that has a random number assigned to each SNP in first column
		cat <(mawk '/^#/' $VCF_FILE) <(paste <(mawk '!/#/' $VCF_FILE | cut -f1-6) $VCF_OUT.mi <(mawk '!/#/' $VCF_FILE | cut -f8- ) )> $VCF_OUT.MostInformativeSNP.temp
		#Use awk (mawk) to parse file and select one snp per contig (one with largest random number)
		cat $VCF_OUT.MostInformativeSNP.temp | mawk 'BEGIN{last_loc = 0} { 
			if ($1 ~/#/) print $0;
			else if ($1 ~!/#/ && last_loc == 0) {last_contig=$0; last_loc=$1; last_qual=$7}
			else if ($1 == last_loc && $7 > last_qual) {last_contig=$0; last_loc=$1; last_qual=$7}
			else if ($1 != last_loc) {print last_contig; last_contig=$0; last_loc=$1; last_qual=$7}
		} END{print last_contig}' | mawk 'NF > 0' > $VCF_OUT.MostInformativeSNP.vcf
		#Remove temp file
		#rm $VCF_OUT.MostInformativeSNP.temp
		mawk '!/#/' $VCF_OUT.MostInformativeSNP.vcf  | wc -l 
		#rm $VCF_OUT.mi
		if [[ $PARALLEL == "TRUE" ]]; then 
			bgzip -@ $NumProc -c $VCF_OUT.MostInformativeSNP.vcf > $VCF_OUT.MostInformativeSNP.vcf.gz
			tabix -p vcf $VCF_OUT.MostInformativeSNP.vcf.gz
		fi	
		
		
	#elif [[ $FILTER_ID == "rmContigs" || $Filter_ID == "99" ]]; then
	elif [[ "$Filter_ID" == "86" ]]; then
		echo; echo `date` "---------------------------FILTER 86: Remove contigs -----------------------------"
		VCF_FILE=${VCF_FILE%.*}
		VCF_FILE_2=${VCF_FILE_2%.*}			
		VCF_OUT=${VCF_FILE%.*}_86rmContigs.vcf
		echo $VCF_FILE $VCF_FILE_2 $VCF_OUT 
		if [[ $PARALLEL == "TRUE" ]]; then 
			printf "${VCF_FILE}\n${VCF_FILE_2}\n" | parallel --no-notice "gunzip -c {}.gz > {}"
			RMcontigs=($(comm -23 ${VCF_FILE_2} ${VCF_FILE} | cut -f1))
			parallel --no-notice -k -j ${NumProc} "grep -v {} $VCF_FILE_2" ::: "${RMcontigs[@]}"  > ${VCF_OUT}
			printf "${VCF_FILE}\n${VCF_FILE_2}\n" | parallel --no-notice rm {}
			bgzip -@ $NumProc -c ${VCF_OUT} > ${VCF_OUT}.gz
			tabix -p vcf ${VCF_OUT}.gz
		else
			#This is a script that identifies which contigs had SNPs that were filtered, then filters those contigs
			#in the next line, the first file is the orig and the second is filtered
			RMcontigs=($(comm -23 $VCF_FILE_2 $VCF_FILE | cut -f1))
			parallel --no-notice -k -j ${NumProc} "grep -v {} $VCF_FILE_2" ::: "${RMcontigs[@]}"  > ${VCF_OUT}
		fi
	fi
}


###################################################################################################################
#specific filter functions
###################################################################################################################

function FILTER_VCFTOOLS(){
	# PARALLEL=$1
	# VCF_FILE=$2
	# Filter=$(echo "$3")
	# VCF_OUT=$4
	# DataName=$5
	# CutoffCode=$6
	# NumProc=$7
	
	# echo "     $Filter"
	# echo "     $FILTER_ID"
	# echo "     $CutoffCode"
	# echo "     $BAM_PATH"
	# echo "     $VCF_FILE"
	# echo "     $REF_FILE"
	# echo "     $PopMap"
	# echo "     $CONFIG_FILE"
	# echo "     $HWE_SCRIPT"
	# echo "     $RADHAP_SCRIPT"
	# echo "     $DataName"
	# echo "     $NumProc"
	# echo "     $PARALLEL"
	# echo "     $COUNTER"
	# echo "     $VCF_OUT"
	
	if [[ $PARALLEL == "FALSE" ]]; then 
		echo "     vcftools --vcf $VCF_FILE $Filter --out $VCF_OUT.recode.vcf 2> /dev/null"
		vcftools --vcf $VCF_FILE $Filter --out $VCF_OUT 2> /dev/null
		echo -n "	Sites remaining:	" && mawk '!/#/' $VCF_OUT.recode.vcf | wc -l
		echo -n "	Contigs remaining:	" && mawk '!/#/' $VCF_OUT.recode.vcf | cut -f1 | uniq | wc -l
	else
		echo "     ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc \"tabix -h -R {} $VCF_FILE | vcftools --vcf - $Filter --stdout 2> /dev/null | tail -n +$NumHeaderLines\" 2> /dev/null | cat $VCF_OUT.header.vcf - | bgzip -@ $NumProc -c > $VCF_OUT.recode.vcf.gz"
		tabix -H $VCF_FILE > $VCF_OUT.header.vcf
		NumHeaderLines=$(tabix -H $VCF_FILE | wc -l)
		NumHeaderLines=$(($NumHeaderLines+1))
		ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcftools --vcf - $Filter --stdout 2> /dev/null | tail -n +$NumHeaderLines" 2> /dev/null | cat $VCF_OUT.header.vcf - | bgzip -@ $NumProc -c > $VCF_OUT.recode.vcf.gz
		tabix -f -p vcf $VCF_OUT.recode.vcf.gz
		rm $VCF_OUT.header.vcf
		echo -n "	Sites remaining:	" && ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.recode.vcf.gz | wc -l " | awk -F: '{a+=$1} END{print a}' 
		echo -n "	Contigs remaining:	" && ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.recode.vcf.gz | cut -f1 | uniq | wc -l " | awk -F: '{a+=$1} END{print a}'
	fi
	
	echo ""
}

function FILTER_VCFFILTER(){
	# PARALLEL=$1
	# VCF_FILE=$2
	# Filter=$(echo "$3")
	# VCF_OUT=$4.vcf
	# DataName=$5
	# CutoffCode=$6
	# NumProc=$7
	VCFFIXUP=$1
	VCF_OUT=$VCF_OUT.vcf
	#Filter=\"$Filter\"

	# echo "     $Filter"
	# echo "     $FILTER_ID"
	# echo "     $CutoffCode"
	# echo "     $BAM_PATH"
	# echo "     $VCF_FILE"
	# echo "     $REF_FILE"
	# echo "     $PopMap"
	# echo "     $CONFIG_FILE"
	# echo "     $HWE_SCRIPT"
	# echo "     $RADHAP_SCRIPT"
	# echo "     $DataName"
	# echo "     $NumProc"
	# echo "     $PARALLEL"
	# echo "     $COUNTER"
	# echo "     $VCF_OUT"
	# echo "     $VCFFIXUP"
	
	if [[ $PARALLEL == "FALSE" ]]; then 
		if [[  $VCFFIXUP == "TRUE" ]]; then echo "     vcffilter -s -f \"$Filter\" $VCF_FILE | vcffixup - > $VCF_OUT"; fi
		if [[  $VCFFIXUP == "FALSE" ]]; then echo "     vcffilter -s -f \"$Filter\" $VCF_FILE > $VCF_OUT"; fi
		if [[  $VCFFIXUP == "TRUE" ]]; then vcffilter -s -f "$Filter" $VCF_FILE | vcffixup - > $VCF_OUT ; fi
		if [[  $VCFFIXUP == "FALSE" ]]; then vcffilter -s -f "$Filter" $VCF_FILE > $VCF_OUT; fi
		echo -n "	Sites remaining:	" && mawk '!/#/' $VCF_OUT | wc -l
		echo -n "	Contigs remaining:	" && mawk '!/#/' $VCF_OUT | cut -f1 | uniq | wc -l
	else
		tabix -H $VCF_FILE > ${VCF_OUT%.*}.header.vcf
		NumHeaderLines=$(tabix -H $VCF_FILE | wc -l)
		NumHeaderLines=$(($NumHeaderLines+2))
		if [[  $VCFFIXUP == "TRUE" ]]; then echo "     ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc \"tabix -h -R {} $VCF_FILE | vcffilter -s -f \"$Filter\" | vcffixup - | tail -n +$NumHeaderLines\" | cat ${VCF_OUT%.*}.header.vcf - | bgzip -@ $NumProc -c > ${VCF_OUT}.gz"; fi
		if [[  $VCFFIXUP == "FALSE" ]]; then echo "     ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc \"tabix -h -R {} $VCF_FILE | vcffilter -s -f \"$Filter\" | tail -n +$NumHeaderLines\" | cat ${VCF_OUT%.*}.header.vcf - | bgzip -@ $NumProc -c > ${VCF_OUT}.gz"; fi
		if [[  $VCFFIXUP == "TRUE" ]]; then ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcffilter -s -f \"$Filter\" | vcffixup - | tail -n +$NumHeaderLines" | cat ${VCF_OUT%.*}.header.vcf - | bgzip -@ $NumProc -c > ${VCF_OUT}.gz; fi
		if [[  $VCFFIXUP == "FALSE" ]]; then ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcffilter -s -f \"$Filter\" | tail -n +$NumHeaderLines" | cat ${VCF_OUT%.*}.header.vcf - | bgzip -@ $NumProc -c > ${VCF_OUT}.gz; fi
		tabix -f -p vcf $VCF_OUT.gz
		#mawk '!/#/' $VCF_OUT.gz | wc -l
		rm ${VCF_OUT%.*}.header.vcf
		echo -n "	Sites remaining:	" && ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.gz | wc -l " | awk -F: '{a+=$1} END{print a}' 
		echo -n "	Contigs remaining:	" && ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.gz | cut -f1 | uniq | wc -l " | awk -F: '{a+=$1} END{print a}'

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
	
	echo ${FILTERS[@]} | sed 's/ /\n/g' > $DataName$CutoffCode.filters.txt
	ls $DataName$CutoffCode.*vcf | parallel --no-notice -k -j $NumProc "grep -c '^dDocent_Contig' {} " | sort -r > $DataName$CutoffCode.SNPcnt.txt
	ls $DataName$CutoffCode.*vcf | parallel --no-notice -k -j $NumProc "grep '^dDocent_Contig' {} | cut -f1 | uniq | wc -l" | sort -r > $DataName$CutoffCode.ContigCnt.txt

	#num missing genotypes for a particular snp
	zgrep -P 'dDocent_Contig_1\t' OpihiSK2014.A.25.10.Fltr18.HWE.recode.vcf.gz | sed '1q;d' | sed 's/\t/\n/g' | grep -c '^\.'
	
}


###################################################################################################################
#help info / manual
###################################################################################################################

NAME="$(basename "$0") v$VERSION  -- a program to filter vcf files with RAD data"

SYNOPSIS="$(basename "$0") [filter settings] [input files] [output file prefix] [parallelization]"

read -d '' DESCRIPTION <<"BLOCK"
fltrVCF is a tool to filter VCF files created by dDocentHPC but can probably work with any work flow. 
		The filters can be run in any order.

        Arguments can be controlled from either the command line or a configuration file.  Filter
        thresholds can only be altered in the config.fltr file. Filters are described  and defined
		in the config.fltr file.

        fltrVCF is parallelized where possible, but only runs on one node or computer. MPI is not
        supported.

        fltrVCF requires minor modification to work with dDocent output.  Both filter_hwe_by_pop_HPC.pl 
		and rad_haplotyper.pl can be obtained from cbirdlab on github and are tested with fltrVCF and work. 
BLOCK

read -d '' OPTIONS <<"BLOCK"
[filter settings]
                -f <arg>        if set, controls filters to be run, in order. Argument should be 2
                                 digit numbers separated by spaces. 
								 -f "01 04 02 86"  or  -f 01\ 04\ 02\ 86 
                                 will specify that filters 01, 04, and 02 will be run in succession.
								 Then, filter 86 will remove the contigs that had SNPs filtered by 02.
                                 filter 86 should only be called after filters that remove SNPs. If 
								 filter 86 is the only filter called, it will compare the newest vcf to
								 the penultimate vcf, and remove contigs based upon the differences.
								 Filters are described in the config files. If -f is not set, the
                                 config file is used to determine the filters and order. If -f is
                                 set, it will override the config file. The results of each filter will
								 be saved in a separate vcf file.[]
                -s <arg>        file with filter settings [config.fltr.ind]. Should be used

		[input files]
                -v <arg>        vcf file to be filtered [${b}/TotalRawSNPs.${c}.vcf]
                -b <arg>        path to mapping directory with *.bam files [../mapping]
				-d <arg>		bed file describing complete ref genome. Required only if -P is set. 
								 [${b}/mapped.${c}.bed]
                -g <arg>        reference genome fasta file [${b}/reference.${c}.fasta]
                -c <arg>        cutoff values used for reference genome, used for file naming, 
				                 optional []
                -p <arg>        popmap file to use for defining population affiliation
                                 [${b}/popmap.${c}]
                -w <arg>        filter_hwe perl script [filter_hwe_by_pop_HPC.pl]
                -r <arg>        rad_haplotyper perl script [rad_haplotyper.pl v1.19]

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
                        -r rad_haplotyper.pl -o ProjectX.A -t 40

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
echo "	rad_haplotyper.pl https://github.com/cbirdlab/rad_haplotyper.git"
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
        CutoffCode=.$OPTARG >&2
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
if [[ -z ${CONFIG_FILE+x} ]]; then 
	CONFIG_FILE=$(ls config.fltr.clean.ind)
	echo "	Settings are being loaded from: '${CONFIG_FILE}' by default"
else 
	CONFIG_FILE=$(ls ${CONFIG_FILE}) 
	echo "	SETTINGS are being loaded from file: '${CONFIG_FILE}'"
fi
if [[ $CONFIG_FILE == "" ]]; then 
	echo "	ERROR :-< 	configuration file is missing" >&2
	exit
fi

if [[ -z ${FILTERS+x} ]]; then 
	FILTERS=($(grep -P '^\t* *fltrVCF\t* *-f\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	echo "	Filters are set to '${FILTERS[@]}'"
else 
	echo "	FILTERS are set to '${FILTERS[@]}'"
fi
if [[ ${FILTERS[0]} == "" ]]; then 
	echo "	ERROR :-< 	filter settings (-f) are missing.  Please modify -f argument in config file or at the command line" >&2
	exit
fi

if [[ -z ${CutoffCode+x} ]]; then 
	CutoffCode=($(grep -P '^\t* *fltrVCF\t* *-c\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $CutoffCode == "" ]]; then
		CutoffCode="" 
	else
		CutoffCode=.$CutoffCode
	fi
fi
echo "	CutoffCode is set to '${CutoffCode}'"

if [[ -z ${BAM_PATH+x} ]]; then 
	BAM_PATH=($(grep -P '^\t* *fltrVCF\t* *-b\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $BAM_PATH == "" ]]; then
		BAM_PATH="../mapping" 
	fi
fi
echo "	BAM_PATH is set to '${BAM_PATH}'"

if [[ -z ${VCF_FILE+x} ]]; then 
	VCF_FILE=($(grep -P '^\t* *fltrVCF\t* *-v\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $VCF_FILE == "" ]]; then
		VCF_FILE=${BAM_PATH}"/TotalRawSNPs"${CutoffCode}".vcf"
	fi
	if [[ ! -f $VCF_FILE ]]; then
		echo "	ERROR:-<	$VCF_FILE does not exist"
		exit
	fi
fi
echo "	VCF_FILE is set to '${VCF_FILE}'"

if [[ -z ${BED_FILE+x} ]]; then 
	BED_FILE=($(grep -P '^\t* *fltrVCF\t* *-d\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ -z ${BED_FILE} ]]; then
		BED_FILE=${BAM_PATH}"/mapped"${CutoffCode}".bed" 
	fi
fi	
echo "	Bed file is set to '${BED_FILE}'"

if [[ -z ${REF_FILE+x} ]]; then 
	REF_FILE=($(grep -P '^\t* *fltrVCF\t* *-g\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $REF_FILE == "" ]]; then
		REF_FILE=${BAM_PATH}"/reference"${CutoffCode}".fasta" 
	fi
fi	
echo "	Reference genome is set to '${REF_FILE}'"

if [[ -z ${PopMap+x} ]]; then 
	PopMap=($(grep -P '^\t* *fltrVCF\t* *-p\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $PopMap == "" ]]; then
		PopMap=${BAM_PATH}"/popmap"${CutoffCode} 
	fi
fi
echo "	PopMap is set to '${PopMap}'"

if [[ -z ${HWE_SCRIPT+x} ]]; then 
	HWE_SCRIPT=($(grep -P '^\t* *fltrVCF\t* *-w\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $HWE_SCRIPT == "" ]]; then
		HWE_SCRIPT="filter_hwe_by_pop_HPC.pl" 
	fi
fi
echo "	HWE_SCRIPT is set to '${HWE_SCRIPT}'"

if [[ -z ${RADHAP_SCRIPT+x} ]]; then 
	RADHAP_SCRIPT=($(grep -P '^\t* *fltrVCF\t* *-r\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $RADHAP_SCRIPT == "" ]]; then
		RADHAP_SCRIPT="rad_haplotyper116HPC.pl" 
	fi
fi
echo "	RADHAP_SCRIPT is set to '${RADHAP_SCRIPT}'"

if [[ -z ${DataName+x} ]]; then 
	DataName=($(grep -P '^\t* *fltrVCF\t* *-o\t* *' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $DataName == "" ]]; then
		echo "	No output file prefix is specified.'" 
	else
		echo "	Output file prefix is set to '${DataName}'"
	fi
else
	echo "	Output file prefix is set to '${DataName}'"
fi


if [[ -z ${NumProc+x} ]]; then 
	NumProc=($(grep -P '^\tfltrVCF -t ' ${CONFIG_FILE} | sed 's/\t* *fltrVCF\t* *-.\t* *//g' | sed 's/\t* *#.*//g')) 
	if [[ $NumProc == "" ]]; then
		NumProc=1 
	fi
fi
echo "	The number of threads is set to '${NumProc}'"

if [[ ${FILTERS[0]} == "" ]]; then 
	echo "	ERROR :-< 	filter settings (-f) are missing.  Please modify -f argument in config file or at the command line" >&2
	exit
fi

if [[ ${PARALLEL} == "TRUE" ]]; then
	if [[ ! -f ${BED_FILE} ]]; then
		echo "	ERROR :-\						Could not find *.bed file. Check the following settings -c -b -d -P."
		echo "									For help, try $(basename "$0") -h"
		exit
	fi
	
	#get contigs with non-overlapping reads into 1 bed file and those with overlapping reads into another
	# cut -f1 ${BED_FILE} | uniq -d > $DataName$CutoffCode.nonoverlapping.contigs
	# cat $DataName$CutoffCode.nonoverlapping.contigs | parallel --no-notice -k -j $NumProc grep {} ${BED_FILE} > $DataName$CutoffCode.nonoverlapping.contigs.bed
	# cut -f1 ${BED_FILE} | uniq -u > $DataName$CutoffCode.overlapping.contigs
	# cat $DataName$CutoffCode.overlapping.contigs | parallel --no-notice -k -j $NumProc grep {} ${BED_FILE} > $DataName$CutoffCode.overlapping.contigs.bed
	
	#shuffle the lines of the bed file
	#cat $(head -n 2 ${BED_FILE) $(tail -n +3 ${BED_FILE} | shuf) 
	
	#split the bed file by the number of processors, making sure that no contig is split between the files
	NumBedLines=$(($(wc -l ${BED_FILE} | sed 's/ .*//g')/$((${NumProc}-1))))
	# RemainderBedLines=$((${NumBedLines}%2))
	# if [[ ${RemainderBedLines} != 0 ]]; then NumBedLines=$((NumBedLines+1)); fi
	split --lines=${NumBedLines} --numeric-suffixes --suffix-length=4 ${BED_FILE} $DataName$CutoffCode.
	ls $DataName$CutoffCode.[0-9][0-9][0-9][0-9] | parallel --no-notice -j $NumProc mv {} {}.bed
	NumBedFiles=$(ls $DataName$CutoffCode.*.bed | wc -l)
	#ensure that no contigs are split between files 
	Indexes=($(seq -f "%04g" 0 $(($NumBedFiles-1))))
	FirstLinePos=($(ls *bed | parallel --no-notice -k head -n 1 {} | cut -f2 )) #takes 2nd col of first lines of each file, except first file
	for i in $(seq 0 $(($NumBedFiles-1)))
	do
		if [[ ${FirstLinePos[$i]} -ge 20 ]]; then   #ceb this interrogates the position; >20 indicates that the pos is not the beginning of the contig
			j=$(($i-1))
			head -n 1 $DataName$CutoffCode.${Indexes[$i]}.bed >> $DataName$CutoffCode.${Indexes[$j]}.bed #copies first line of bed.index to last line of bed.index-1
			sed -i '1d' $DataName$CutoffCode.${Indexes[$i]}.bed #removes first line of bed.index
		fi
	done
	
	echo "	-P set. All filters will be coerced to run in parallel."
	if [[ "${VCF_FILE##*.}" == "vcf" ]]; then
		if [[ ! -f "${VCF_FILE}.gz" ]]; then
			echo ""; echo "	" `date` " Using bgzip and tabix to compress and index VCF for parallelization"
			bgzip -@ $NumProc -c ${VCF_FILE} > ${VCF_FILE}.gz
			VCF_FILE=${VCF_FILE}.gz 
			tabix -p vcf ${VCF_FILE}
		else
			echo "	Changing -v to existing file ${VCF_FILE}.gz for parallelization.  If you didn't want this to happen, turn off -P or change -v"
			VCF_FILE=${VCF_FILE}.gz
		fi
	elif [[ "${VCF_FILE##*.}" == "gz" ]]; then
		if [[ ! -f "${VCF_FILE}" ]]; then
			echo "	ERROR:-(   this should be impossible to trigger"
		else
			echo "	${VCF_FILE} will be used for parallelization."
		fi
		if [[ ! -f "${VCF_FILE}.tbi" ]]; then
			echo "	${VCF_FILE}.tbi not found. Using tabix to index the *vcf.gz file"
			tabix -p vcf ${VCF_FILE}
		else
			echo "	${VCF_FILE}.tbi was successfully located."
		fi
	fi
else
	echo "	-P not set by user. Only filters that natively support parallelization will be run in parallel."
	PARALLEL=FALSE
	if [[ "${VCF_FILE##*.}" == "gz" ]]; then
		if [[ ! -f "${VCF_FILE%.*}" ]]; then
			echo ""; echo "	" `date` " Using gunzip to decompress $VCF_FILE"
			gunzip ${VCF_FILE} > ${VCF_FILE%.*}
			VCF_FILE=${VCF_FILE%.*} 
		else
			echo "	Changing -v to existing file ${VCF_FILE%.*} for serial processing.  If you didn't want this to happen, turn on -P or change -v"
			VCF_FILE=${VCF_FILE%.*}
		fi
	elif [[ "${VCF_FILE##*.}" == "vcf" ]]; then
		if [[ ! -f "${VCF_FILE}" ]]; then
			echo "	ERROR:-(   this should be impossible to trigger"
		else
			echo "	${VCF_FILE} will be used for serial processing."
		fi
	fi
fi

echo ""

###################################################################################################################
#Run script
###################################################################################################################

MAIN #FILTERS $CutoffCode $BAM_PATH $VCF_FILE $REF_FILE $PopMap $CONFIG_FILE $HWE_SCRIPT $RADHAP_SCRIPT $DataName $NumProc $PARALLEL

#cleanup files
if [[ ${PARALLEL} == "TRUE" ]]; then
	ls $DataName$CutoffCode.[0-9][0-9][0-9][0-9].bed | parallel --no-notice -j $NumProc "rm {}"
fi
echo ""; echo `date` " --------------------------- Filtering complete! ---------------------------"; echo "" 