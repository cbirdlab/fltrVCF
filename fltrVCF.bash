#!/bin/bash

#Improvements yet to be made
#	Shuffle the bed file to improve performance
#	Finish parallel option for rad_haplotyper
#	Enable option to make haplotype files with radhaplotyper
#	figure out how to combine haplotype files made by radhaplotyper in parallel mode

VERSION=4.3
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
	
	if [[ $PARALLEL == "FALSE" ]]; then
#		VCF_FILE=$(ls -t ${DataName}*vcf | head -n 1)
		VCF_FILE_2=$(ls -t ${DataName}*vcf | sed -n 2p)
	else
#		VCF_FILE=$(ls -t ${DataName}*vcf.gz | head -n 1)
		VCF_FILE_2=$(ls -t ${DataName}*vcf.gz | sed -n 2p)
	fi	
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
	if [[ $FILTER_ID == "30" ]]; then
		echo; echo `date` "---------------------------FILTER30: Remove Positions From Consideration -----------------------------"
		COUNTER30=$((1+COUNTER30))
		#get settings from config file
		THRESHOLD=$(grep -P '^\t* *30\t* *custom\t* *bash\t* *..*\t* *#Keep sites after this position' ${CONFIG_FILE} | sed 's/\t* *30\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLD" ]]; then THRESHOLD=0; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		THRESHOLDb=$(grep -P '^\t* *30\t* *custom\t* *bash\t* *..*\t* *#Keep sites before this position' ${CONFIG_FILE} | sed 's/\t* *30\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLDb" ]]; then THRESHOLDb=5000; fi
		THRESHOLDb=$(PARSE_THRESHOLDS $THRESHOLDb) 
		if [[ $PARALLEL == "FALSE" ]]; then
			grep '^dDocent' $VCF_FILE | grep -w "TYPE=snp" | cut -f1-4 > $VCF_OUT.prefltr.snp &
			grep '^dDocent' $VCF_FILE | grep -w "TYPE=ins" | cut -f1-4 > $VCF_OUT.prefltr.ins &
			grep '^dDocent' $VCF_FILE | grep -w "TYPE=del" | cut -f1-4 > $VCF_OUT.prefltr.del &
			vcftools --vcf ${VCF_FILE} --site-mean-depth --out $VCF_OUT.prefltr 2> /dev/null &
			wait
			tail -n+2 $VCF_OUT.prefltr.ldepth.mean | cut -f1-3 | grep -v 'nan' > $VCF_OUT.ldepth.mean.prefltr30
			Rscript plotFltr30.R $VCF_OUT.ldepth.mean.prefltr30 $VCF_OUT.prefltr.snp $VCF_OUT.prefltr.ins $VCF_OUT.prefltr.del $THRESHOLD $THRESHOLDb $VCF_OUT.prefltr.plots.pdf
			echo " Plots output to $VCF_OUT.prefltr.plots.pdf"

			cat <(grep -Pv '^dDocent_Contig_[1-9][0-9]*' ${VCF_FILE}) \
				<(grep -P '^dDocent_Contig_[1-9][0-9]*' ${VCF_FILE} | awk -v bp=$THRESHOLD '$2 > bp {print ;}' | awk -v bp=$THRESHOLDb '$2 < bp {print ;}' ) \
				> $VCF_OUT.vcf
			
			grep '^dDocent' $VCF_OUT.vcf | grep -w "TYPE=snp" | cut -f1-4 > $VCF_OUT.postfltr.snp &
			grep '^dDocent' $VCF_OUT.vcf | grep -w "TYPE=ins" | cut -f1-4 > $VCF_OUT.postfltr.ins &
			grep '^dDocent' $VCF_OUT.vcf | grep -w "TYPE=del" | cut -f1-4 > $VCF_OUT.postfltr.del &
			vcftools --vcf $VCF_OUT.vcf --site-mean-depth --out $VCF_OUT.postfltr 2> /dev/null &
			wait
			tail -n+2 $VCF_OUT.postfltr.ldepth.mean | cut -f1-3 | grep -v 'nan' > $VCF_OUT.ldepth.mean.postfltr30
			Rscript plotFltr30.R $VCF_OUT.ldepth.mean.postfltr30 $VCF_OUT.postfltr.snp $VCF_OUT.postfltr.ins $VCF_OUT.postfltr.del $THRESHOLD $THRESHOLDb $VCF_OUT.postfltr.plots.pdf
			
		else
			zgrep '^dDocent' $VCF_FILE | grep -w "TYPE=snp" | cut -f1-4 > $VCF_OUT.prefltr.snp &
			zgrep '^dDocent' $VCF_FILE | grep -w "TYPE=ins" | cut -f1-4 > $VCF_OUT.prefltr.ins &
			zgrep '^dDocent' $VCF_FILE | grep -w "TYPE=del" | cut -f1-4 > $VCF_OUT.prefltr.del &
			vcftools --vcf ${VCF_FILE} --site-mean-depth --out $VCF_OUT.prefltr 2> /dev/null &
			#ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -h -R {} $VCF_FILE | vcftools --vcf - $Filter --stdout 2> /dev/null | tail -n +$NumHeaderLines" 2> /dev/null | cat $VCF_OUT.header.vcf - | bgzip -@ $NumProc -c > $VCF_OUT.recode.vcf.gz
			wait
			tail -n+2 $VCF_OUT.prefltr.ldepth.mean | cut -f1-3 | grep -v 'nan' > $VCF_OUT.ldepth.mean.prefltr.fltr30
			Rscript plotFltr30.R $VCF_OUT.ldepth.mean.prefltr.fltr30 $VCF_OUT.prefltr.snp $VCF_OUT.prefltr.ins $VCF_OUT.prefltr.del $THRESHOLD $THRESHOLDb $VCF_OUT.prefltr.plots.pdf
			echo " Plots output to $VCF_OUT.prefltr.plots.pdf"
			cat <(zgrep -Pv '^dDocent_Contig_[1-9][0-9]*' ${VCF_FILE}) \
				<(zgrep -P '^dDocent_Contig_[1-9][0-9]*' ${VCF_FILE} | awk -v bp=$THRESHOLD '$2 > bp {print ;}' | awk -v bp=$THRESHOLDb '$2 < bp {print ;}' ) \
				| bgzip -@ $NumProc -c > $VCF_OUT.vcf.gz
			tabix -f -p vcf $VCF_OUT.vcf.gz

			zgrep '^dDocent' $VCF_OUT.vcf.gz | grep -w "TYPE=snp" | cut -f1-4 > $VCF_OUT.postfltr.snp &
			zgrep '^dDocent' $VCF_OUT.vcf.gz | grep -w "TYPE=ins" | cut -f1-4 > $VCF_OUT.postfltr.ins &
			zgrep '^dDocent' $VCF_OUT.vcf.gz | grep -w "TYPE=del" | cut -f1-4 > $VCF_OUT.postfltr.del &
			vcftools --vcf $VCF_OUT.vcf.gz --site-mean-depth --out $VCF_OUT.postfltr 2> /dev/null &
			wait
			tail -n+2 $VCF_OUT.postfltr.ldepth.mean | cut -f1-3 | grep -v 'nan' > $VCF_OUT.ldepth.mean.postfltr30
			Rscript plotFltr30.R $VCF_OUT.ldepth.mean.postfltr30 $VCF_OUT.postfltr.snp $VCF_OUT.postfltr.ins $VCF_OUT.postfltr.del $THRESHOLD $THRESHOLDb $VCF_OUT.postfltr.plots.pdf
		fi
		echo " Plots output to $VCF_OUT.postfltr.plots.pdf"


	elif [[ $FILTER_ID == "31" ]]; then
		echo; echo `date` "---------------------------FILTER31: Remove Contigs With Fewer BP -----------------------------"
		COUNTER31=$((1+COUNTER31))
		#get settings from config file
		THRESHOLD=$(grep -P '^\t* *31\t* *custom\t* *bash\t* *..*\t* *#Remove contigs with fewer BP' ${CONFIG_FILE} | sed 's/\t* *31\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLD" ]]; then THRESHOLD=0; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 

		if [[ $PARALLEL == "FALSE" ]]; then
			#calculate the mean depth by site with VCFtools
			vcftools --vcf ${VCF_FILE} --site-mean-depth --out $VCF_OUT 2> /dev/null
			tail -n+2 $VCF_OUT.ldepth.mean | cut -f1-3 | grep -v 'nan' | awk '{x[$1] += $3; N[$1]++} END{for (i in x) print i, N[i], x[i]/N[i]}' | tr -s " " "\t" | sort -nk3 > $VCF_OUT.ldepth.mean.contigs
			Rscript plotFltr31.R $VCF_OUT.ldepth.mean.contigs $THRESHOLD $VCF_OUT.ldepth.mean.contigs.plots.pdf
			awk -v BP=$THRESHOLD '$2 <= BP {print $1;}' $VCF_OUT.ldepth.mean.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			grep -vf $VCF_OUT.remove.contigs $VCF_FILE > $VCF_OUT.vcf
		else
			#calculate the mean depth by site with VCFtools
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			tail -n+2 $VCF_OUT.ldepth.mean | cut -f1-3 | grep -v 'nan' | awk '{x[$1] += $3; N[$1]++} END{for (i in x) print i, N[i], x[i]/N[i]}' | tr -s " " "\t" | sort -nk3 > $VCF_OUT.ldepth.mean.contigs
			Rscript plotFltr31.R $VCF_OUT.ldepth.mean.contigs $THRESHOLD $VCF_OUT.ldepth.mean.contigs.plots.pdf
			awk -v BP=$THRESHOLD '$2 <= BP {print $1;}' $VCF_OUT.ldepth.mean.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			zgrep -vf $VCF_OUT.remove.contigs $VCF_FILE | bgzip -@ $NumProc -c > $VCF_OUT.vcf.gz
			tabix -f -p vcf $VCF_OUT.vcf.gz
		fi
		echo " The following contigs have been filtered:"
		sed -e 's/^\^/\t/' -e 's/\t$//' $VCF_OUT.remove.contigs | sort -g | paste -d "" - - - - -
		echo ""; echo " Plots output to $VCF_OUT.ldepth.mean.contigs.plots.pdf"

	elif [[ $FILTER_ID == "32" ]]; then
		echo; echo `date` "---------------------------FILTER32: Keep Contigs w Fewer Heterozygotes -----------------------------"
		COUNTER32=$((1+COUNTER32))
		#get settings from config file
		THRESHOLD=$(grep -P '^\t* *32\t* *custom\t* *bash\t* *..*\t* *#Keep contigs with lesser porportion of heterozygotes' ${CONFIG_FILE} | sed 's/\t* *32\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLD" ]]; then THRESHOLD=1; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 

		if [[ $PARALLEL == "FALSE" ]]; then
			#get heterozygosity counts, then get mean rate of heterozygosity per variable position per contig, accounting for missing data
			grep '^dDocent' ${VCF_FILE} | awk 'BEGIN{print "Contig\tNumHet\tNumHomoRef\tNumHomoAlt\tNumMissing\tNumInd"}{print $1 "\t" gsub(/0\/1:/,"0/1:") "\t" gsub(/0\/0:/,"0/0:") "\t" gsub(/1\/1:/,"1/1:") "\t" gsub(/\.\/\.:/,"./.:") "\t" gsub(/.\/.:/,"") } ' | \
				tail -n+2 | awk '$3 + $5 != $6 && $4 + $5 != $6 {print ;}' | \
				awk '{x[$1] += $2; y[$1] += $6 -=$5; N[$1]++} END{for (i in x) print i, N[i], x[i]/N[i], y[i]/N[i], x[i]/y[i] }' | tr -s " " "\t" | sort -nrk5 > $VCF_OUT.hetero.contigs
			awk -v HET=$THRESHOLD '$5 >= HET {print $1;}' $VCF_OUT.hetero.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			Rscript plotFltr32.R $VCF_OUT.hetero.contigs $THRESHOLD $VCF_OUT.hetero.contigs.plots.pdf
			grep -vf $VCF_OUT.remove.contigs $VCF_FILE > $VCF_OUT.vcf
		else
			#get heterozygosity counts, then get mean rate of heterozygosity per variable position per contig, accounting for missing data
			zgrep '^dDocent' ${VCF_FILE} | awk 'BEGIN{print "Contig\tNumHet\tNumHomoRef\tNumHomoAlt\tNumMissing\tNumInd"}{print $1 "\t" gsub(/0\/1:/,"0/1:") "\t" gsub(/0\/0:/,"0/0:") "\t" gsub(/1\/1:/,"1/1:") "\t" gsub(/\.\/\.:/,"./.:") "\t" gsub(/.\/.:/,"") } ' | \
				tail -n+2 | awk '$3 + $5 != $6 && $4 + $5 != $6 {print ;}' | \
				awk '{x[$1] += $2; y[$1] += $6 -=$5; N[$1]++} END{for (i in x) print i, N[i], x[i]/N[i], y[i]/N[i], x[i]/y[i] }' | tr -s " " "\t" | sort -nrk5 > $VCF_OUT.hetero.contigs
			awk -v HET=$THRESHOLD '$5 >= HET {print $1;}' $VCF_OUT.hetero.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			Rscript plotFltr32.R $VCF_OUT.hetero.contigs $THRESHOLD $VCF_OUT.hetero.contigs.plots.pdf
			zgrep -vf $VCF_OUT.remove.contigs $VCF_FILE | bgzip -@ $NumProc -c > $VCF_OUT.vcf.gz
			tabix -f -p vcf $VCF_OUT.vcf.gz
		fi
		echo " The following contigs have been filtered:"
		sed -e 's/^\^/\t/' -e 's/\t$//' $VCF_OUT.remove.contigs | sort -n | paste -d "" - - - - -
		echo ""; echo " Plots output to $VCF_OUT.hetero.contigs.plots.pdf"	
		
	elif [[ $FILTER_ID == "041" ]]; then
		echo; echo `date` "---------------------------FILTER041: Remove Contigs With Extreme DP -----------------------------"
		COUNTER041=$((1+COUNTER041))
		#get settings from config file
THRESHOLD=$(grep -P '^\t* *041\t* *custom\t* *bash\t* *..*\t* *#Remove contigs with lower mean of mean depth across sites' ${CONFIG_FILE} | sed 's/\t* *041\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLD" ]]; then THRESHOLD=0.01; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		THRESHOLDb=$(grep -P '^\t* *041\t* *custom\t* *bash\t* *..*\t* *#Remove contigs with higher mean of mean depth across sites' ${CONFIG_FILE} | sed 's/\t* *041\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLDb" ]]; then THRESHOLDb=0.99; fi
		THRESHOLDb=$(PARSE_THRESHOLDS $THRESHOLDb) 
		THRESHOLDc=$(grep -P '^\t* *041\t* *custom\t* *bash\t* *..*\t* *#Remove contigs with lower CV of mean depth across sites' ${CONFIG_FILE} | sed 's/\t* *041\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLDc" ]]; then THRESHOLDc=0; fi
		THRESHOLDc=$(PARSE_THRESHOLDS $THRESHOLDc) 
		THRESHOLDd=$(grep -P '^\t* *041\t* *custom\t* *bash\t* *..*\t* *#Remove contigs with higher CV of mean depth across sites' ${CONFIG_FILE} | sed 's/\t* *041\t* *custom\t* *bash\t* *//g' | sed 's/\t* *#.*//g' ) 
		if [[ -z "$THRESHOLDd" ]]; then THRESHOLDd=0.99; fi
		THRESHOLDd=$(PARSE_THRESHOLDS $THRESHOLDd) 

		if [[ $PARALLEL == "FALSE" ]]; then
			#calculate the mean depth by site with VCFtools
			vcftools --vcf ${VCF_FILE} --site-mean-depth --out $VCF_OUT 2> /dev/null
			#this awk command calculates the number of bp considered in the contig, the mean of mean cvg, the var or mean cvg, the mean of CV for the mean cvg, and the var of the CV for the mean cvg
			tail -n+2 $VCF_OUT.ldepth.mean | grep -Pv '\t-nan' | awk '{sumAVG[$1] += $3; ssAVG[$1] += $3 ** 2; sumCV[$1] += $4 ** 0.5 / $3; ssCV[$1] += ($4 ** 0.5 / $3)**2; N[$1]++} END{for (i in sumAVG) print i, N[i], sumAVG[i]/N[i], (ssAVG[i]/N[i] - (sumAVG[i]/N[i])**2)**0.5 / sumAVG[i]/N[i], sumCV[i]/N[i], (ssCV[i]/N[i] - (sumCV[i]/N[i])**2)**0.5 / (sumCV[i]/N[i]) }' | tr -s " " "\t" | sort -nk3 > $VCF_OUT.ldepth.mean.contigs
			Rscript plotFltr041.R $VCF_OUT.ldepth.mean.contigs $THRESHOLD $THRESHOLDb $THRESHOLDc $THRESHOLDd $VCF_OUT.ldepth.mean.contigs.plots.pdf
			THRESHOLD=$(cut -f3 $VCF_OUT.ldepth.mean.contigs | awk -v PCT=$THRESHOLD '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			THRESHOLDb=$(cut -f3 $VCF_OUT.ldepth.mean.contigs | awk -v PCT=$THRESHOLDb '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			THRESHOLDc=$(cut -f5 $VCF_OUT.ldepth.mean.contigs | sort -n | awk -v PCT=$THRESHOLDc '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			if [[ -z "$THRESHOLDc" ]]; then THRESHOLDc=0; fi
			THRESHOLDd=$(cut -f5 $VCF_OUT.ldepth.mean.contigs | sort -n | awk -v PCT=$THRESHOLDd '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			#this was line before adding cv filter, delete after confirmation of success
			#awk -v CVGb=$THRESHOLDb -v CVG=$THRESHOLD '$3 > CVGb || $3 < CVG {print $1;}' $VCF_OUT.ldepth.mean.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			awk -v CVGd=$THRESHOLDd -v CVGc=$THRESHOLDc -v CVGb=$THRESHOLDb -v CVG=$THRESHOLD '$3 > CVGb || $3 < CVG || $5 > CVGd || $5 < CVGc {print $1;}' $VCF_OUT.ldepth.mean.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			grep -vf $VCF_OUT.remove.contigs $VCF_FILE > $VCF_OUT.vcf
		else
			#calculate the mean depth by site with VCFtools
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			tail -n+2 $VCF_OUT.ldepth.mean | cut -f1-3 | grep -v 'nan' | awk '{meancvg[$1] += $3; N[$1]++} END{for (i in meancvg) print i, N[i], meancvg[i]/N[i]}' | tr -s " " "\t" | sort -nk3 > $VCF_OUT.ldepth.mean.contigs
			Rscript plotFltr041.R $VCF_OUT.ldepth.mean.contigs $THRESHOLD $THRESHOLDb $THRESHOLDc $THRESHOLDd $VCF_OUT.ldepth.mean.contigs.plots.pdf
			THRESHOLD=$(cut -f3 $VCF_OUT.ldepth.mean.contigs | awk -v PCT=$THRESHOLD '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			THRESHOLDb=$(cut -f3 $VCF_OUT.ldepth.mean.contigs | awk -v PCT=$THRESHOLDb '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			THRESHOLDc=$(cut -f5 $VCF_OUT.ldepth.mean.contigs | sort -n | awk -v PCT=$THRESHOLDc '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			THRESHOLDd=$(cut -f5 $VCF_OUT.ldepth.mean.contigs | sort -n | awk -v PCT=$THRESHOLDd '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
			awk -v CVGd=$THRESHOLDd -v CVGc=$THRESHOLDc -v CVGb=$THRESHOLDb -v CVG=$THRESHOLD '$3 > CVGb || $3 < CVG || $5 > CVGd || $5 < CVGc {print $1;}' $VCF_OUT.ldepth.mean.contigs | sed -e 's/^/\^/' -e 's/$/\t/' > $VCF_OUT.remove.contigs
			zgrep -vf $VCF_OUT.remove.contigs $VCF_FILE | bgzip -@ $NumProc -c > $VCF_OUT.vcf.gz
			tabix -f -p vcf $VCF_OUT.vcf.gz
		fi

		echo " The following contigs have been filtered:"
		sed -e 's/^\^/\t/' -e 's/\t$//' $VCF_OUT.remove.contigs | sort -g | paste -d "" - - - - -
		echo ""; echo " Plots output to $VCF_OUT.ldepth.mean.contigs.pdf"
		
	elif [[ $FILTER_ID == "01" ]]; then
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
		if [[ $PARALLEL == "FALSE" ]]; then 
			grep '^dDocent_Contig' ${VCF_FILE} | cut -f6 | sort -rg > ${VCF_OUT}.QUAL.before.dat
		else
			zgrep '^dDocent_Contig' ${VCF_FILE} | cut -f6 | sort -rg > ${VCF_OUT}.QUAL.before.dat
		fi
		cp ${VCF_OUT}.QUAL.before.dat QUALbefore
		
gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale
	set xrange [0:*] 
	unset label
	set title "Histogram of QUAL before FILTER03"
	set ylabel "Number of Sites"
	set xlabel "QUAL ~ Phred-Scaled Probality of Zero Alternate Alleles"
	xmax=500
	set xrange [0:xmax]
	binwidth=xmax/50
	bin(x,width)=width*floor(x/width) + binwidth/2.0
	#set xtics 5
	plot 'QUALbefore' using (bin($1,binwidth)):(1.0) smooth freq with boxes
	#plot 'QUALbefore' smooth freq with boxes
	pause -1
EOF

gnuplot << \EOF 
	set terminal dumb size 120, 30
	set autoscale 
	unset label
	set title "Scatter plot of QUAL per site."
	set ylabel "QUAL ~ Phred-Scaled Probality of Zero Alternate Alleles"
	set xlabel "Site"
	#xmax="`cut -f1 QUALbefore | tail -1`"
	#xmax=xmax+1
	#set xrange [0:xmax]
	set yrange [0:500]
	plot 'QUALbefore' pt "*" 
	pause -1
EOF
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "04" ]]; then
		echo; echo `date` "---------------------------FILTER04: Remove Sites With Mean Depth of Coverage < min-meanDP -----------------------------"
		THRESHOLD=($(grep -P '^\t* *04\t* *vcftools\t* *--min-meanDP' ${CONFIG_FILE} | sed 's/\t* *04\t* *vcftools\t* *--min-meanDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=2; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter2="--site-mean-depth"
		Filter="--min-meanDP ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		
		#calculate the mean depth by site with VCFtools
		if [[ $PARALLEL == "FALSE" ]]; then 
			vcftools --vcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			cut -f3 $VCF_OUT.ldepth.mean > $VCF_OUT.site.depth.mean
		else
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			cut -f3 $VCF_OUT.ldepth.mean > $VCF_OUT.site.depth.mean
		fi
		tail -n +2 $VCF_OUT.site.depth.mean | sort -rg > sitedepthmeanbefore
		cut -f3-4 $VCF_OUT.ldepth.mean | tail -n +2 > meandepthVSvariance
		
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
set xrange [0:*] 
unset label
set title "Histogram of mean depth per site before FILTER04"
set ylabel "Number of Sites"
set xlabel "Mean Depth"
binwidth=1
bin(x,width)=width*floor(x/width) + binwidth/2.0
#set xtics 5
plot 'sitedepthmeanbefore' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of mean depth per site before FILTER04."
set ylabel "Mean Depth"
set xlabel "Site"
plot 'sitedepthmeanbefore' pt "*" 
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of depth variance vs mean depth per site before FILTER04."
set ylabel "Variance in Depth"
set xlabel "Mean Depth"
plot 'meandepthVSvariance' pt "*"
pause -1
EOF
		
		FILTER_VCFTOOLS #$PARALLEL $VCF_FILE "${Filter}" $VCF_OUT $DataName $CutoffCode $NumProc 

	elif [[ $FILTER_ID == "05" ]]; then
		echo; echo `date` "---------------------------FILTER05: Remove sites called in <X proportion of individuals -----------------------------"
		THRESHOLD=($(grep -P '^\t* *05\t* *vcftools\t* *--max-missing' ${CONFIG_FILE} | sed 's/\t* *05\t* *vcftools\t* *--max-missing\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.5; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter2="--missing-site"
		Filter="--max-missing ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID

		if [[ $PARALLEL == "FALSE" ]]; then 
			vcftools --vcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
		else
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
		fi
		cut -f6 $VCF_OUT.lmiss | tail -n +2 | sort -rg > sitemissingness
		
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of the proportion of genotypes that are missing"
set ylabel "Number of Sites"
set xlabel "Proportion Missing"
set xrange [0:1]
binwidth=0.05
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'sitemissingness' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of the proportion of genotypes that are missing."
set ylabel "Proportion Missing Genotypes"
set xlabel "Site"
plot 'sitemissingness' pt "*" 
pause -1
EOF


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
		
		vcf-query ${VCF_FILE}  -f '%INFO/AB\n' > $VCF_OUT.siteAB.dat
		sort -rg $VCF_OUT.siteAB.dat > siteAB.dat
		
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of Allele Balance"
set ylabel "Number of Sites"
set xlabel "Allele Balance ~ proportion of alleles in heterozygotes that are reference allelic state"
set xrange [0:1]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'siteAB.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Allele Balance."
set ylabel "Allele Balance"
set xlabel "Site"
plot 'siteAB.dat' pt "*" 
pause -1
EOF

		
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

		# if [[ $PARALLEL == "FALSE" ]]; then 
			# vcftools --vcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			# #Now let’s take VCFtools output and cut it to only the depth scores
			# cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			# #Now let’s calculate the average depth by dividing the above file by the number of individuals
			# NumInd=$(awk '{if ($1 == "#CHROM"){print NF-9; exit}}' $VCF_FILE )
		# else
			# vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
			# #Now let’s take VCFtools output and cut it to only the depth scores
			# cut -f3 $VCF_OUT.ldepth > $VCF_OUT.site.depth
			# #Now let’s calculate the average depth by dividing the above file by the number of individuals
			# NumInd=$(gunzip -c $VCF_FILE | awk '{if ($1 == "#CHROM"){print NF-9; exit}}' )
		# fi
		# mawk '!/D/' $VCF_OUT.site.depth | mawk -v x=$NumInd '{print $1/x}' > $VCF_OUT.meandepthpersiteafter
		# cp $VCF_OUT.meandepthpersiteafter meandepthpersiteafter
		
# gnuplot << \EOF 
# set terminal dumb size 120, 30
# set autoscale
# set xrange [0:*] 
# unset label
# set title "Histogram of mean depth per site after filter"
# set ylabel "Number of Occurrences"
# set xlabel "Mean Depth"
# binwidth=1
# bin(x,width)=width*floor(x/width) + binwidth/2.0
# #set xtics 5
# plot 'meandepthpersiteafter' using (bin($1,binwidth)):(1.0) smooth freq with boxes
# pause -1
# EOF
		
	elif [[ $FILTER_ID == "14" ]]; then
		echo; echo `date` "---------------------------FILTER14: If individual's genotype has DP < X, convert to missing data -----------------------------"
		THRESHOLD=($(grep -P '^\t* *14\t* *vcftools\t* *--minDP' ${CONFIG_FILE} | sed 's/\t* *14\t* *vcftools\t* *--minDP\t* *//g' | sed 's/\t* *#.*//g' )) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=3; fi
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		Filter2="--geno-depth"
		Filter="--minDP ${THRESHOLD} --recode --recode-INFO-all"
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		
				#calculate the mean depth by site with VCFtools
		if [[ $PARALLEL == "FALSE" ]]; then 
			vcftools --vcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
		else
			vcftools --gzvcf ${VCF_FILE} $Filter2 --out $VCF_OUT 2> /dev/null
		fi
		cut -f3- $VCF_OUT.gdepth | tail -n +2 > $VCF_OUT.genotype.depth
		NumCols=$(head -1 $VCF_OUT.genotype.depth | tr "\t" "\n" | wc -l)
		seq 1 $NumCols | parallel --no-notice -j $NumProc "cut -f {} $VCF_OUT.genotype.depth " | sort -rg | sed 's/-1/0/g' > genotypedepthbefore
		
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale
unset label
set title "Histogram of genotype depth before FILTER14"
set ylabel "Number of Genotypes"
set xlabel "Depth"
xmax=100
set xrange [0:xmax]
binwidth=xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'genotypedepthbefore' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of mean depth per site before FILTER14."
set ylabel "Depth"
set xlabel "Genotype"
set yrange [0:100]
plot 'genotypedepthbefore' pt "*" 
pause -1
EOF


		
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
		grep -P '^#CHROM\tPOS' $VCF_OUT.recode.vcf > $VCF_OUT.header.line
		# while read i; do
			# sed -i "s/\t$i//g" $VCF_OUT.header.line
		# done < $VCF_OUT.lowDP-2.indv
		for i in $VCF_OUT.lowDP-2.indv; do
			sed -i "s/\t$i//g" $VCF_OUT.header.line
		done
		#replace header line in vcf
		sed -i "s/^#CHROM\tPOS\t.*$/$(cat $VCF_OUT.header.line)/" $VCF_OUT.recode.vcf
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
		
		#while re
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
		
		
	elif [[ $FILTER_ID == "181" ]]; then
		echo; echo `date` "---------------------------FILTER181: Remove sites & contigs in & not in HWE p<X) -----------------------------"
		THRESHOLD=($(grep -P '^\t* *18\t* *filter_hwe_by_pop_HPC' ${CONFIG_FILE} | sed 's/\t* *18\t* *filter_hwe_by_pop_HPC\t* *//g' | sed 's/\t* *#.*//g' )) 
		THRESHOLD=$(PARSE_THRESHOLDS $THRESHOLD) 
		if [[ -z "${THRESHOLD}" ]]; then ${THRESHOLD}=0.001; fi
		#VCF_OUT=$DataName$CutoffCode.Fltr$FILTER_ID
		if [[ $PARALLEL == "TRUE" ]]; then 
			gunzip -c $VCF_FILE > ${VCF_FILE%.*}
			VCF_FILE=${VCF_FILE%.*}
		fi
		perl $HWE_SCRIPT -v $VCF_FILE -p $PopMap -h ${THRESHOLD} -d $DataName -co $CutoffCode -o $VCF_OUT.HWE
		cat <(grep -v '^dDocent_Contig' ${VCF_FILE}) <(comm -23 <(grep '^dDocent_Contig' ${VCF_FILE} | sort -g) <(grep '^dDocent_Contig' $VCF_OUT.HWE.recode.vcf | sort -g) | sort -V ) > $VCF_OUT.SITES.NOT.IN.HWE.vcf
		cat $VCF_OUT.SITES.NOT.IN.HWE.vcf | grep '^dDocent_Contig' | cut -f1 | uniq > contigs_failing_HWE.txt
		grep -f "contigs_failing_HWE.txt" ${VCF_FILE} > $VCF_OUT.CONTIGS.NOT.IN.HWE.vcf
		grep -v -f "contigs_failing_HWE.txt" ${VCF_FILE} > $VCF_OUT.CONTIGS.IN.HWE.vcf
		mawk '!/#/' $VCF_OUT.CONTIGS.IN.HWE.vcf | wc -l
		if [[ $PARALLEL == "TRUE" ]]; then 
			bgzip -@ $NumProc -c $VCF_OUT.CONTIGS.IN.HWE.vcf > $VCF_OUT.SITES.NOT.IN.HWE.vcf.gz
			tabix -p vcf $VCF_OUT.CONTIGS.IN.HWE.vcf.gz
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
					perl $RADHAP_SCRIPT -v ${Indexing}.vcf -x 1 -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}$CutoffCode -o ${VCF_OUT}.Fltr${FILTER_ID}.Haplotyped.vcf -g ${VCF_OUT}.Fltr${FILTER_ID}.${PopMap##*/}.haps.genepop -a ${VCF_OUT}.Fltr${FILTER_ID}.${PopMap##*/}.haps.ima

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
			
			perl $RADHAP_SCRIPT -v $VCF_FILE.1.vcf -x ${NumProc} -e -d ${THRESHOLDa} -mp ${THRESHOLDb} -u ${THRESHOLDc} -ml ${THRESHOLDd} -h ${THRESHOLDe} -z ${THRESHOLDf} -m ${THRESHOLDg} -r ${REF_FILE} -bp ${BAM_PATH} -p ${PopMap}$CutoffCode -o $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf -g $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.genepop -a $VCF_OUT.Fltr$FILTER_ID.${PopMap##*/}.haps.ima
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
		echo ""
		Explanations=("paralog" "Missing" "haplotypes" "Coverage" "SNPs" "Complex")
		parallel --gnu --null "grep {} $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f1 > $VCF_OUT.Fltr$FILTER_ID.{}" ::: "${Explanations[@]}"
		parallel --gnu --null "echo -n {}' removed, ' && wc -l $VCF_OUT.Fltr$FILTER_ID.{}" ::: "${Explanations[@]}"
		
#histogram of # SNPs
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f2 | sort -rg > sites.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of number of sites per contig."
set ylabel "Number of Contigs"
set xlabel "Number of Sites"
set yrange [0:*]
xmax="`head -1 sites.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'sites.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of # Sites per Contig."
set ylabel "# Sites"
set xlabel "Contigs"
plot 'sites.dat' pt "*" 
pause -1
EOF

#histogram of # haplotypes
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f3 | sort -rg > haplotypes.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of number of haplotypes per contig."
set ylabel "Number of Contigs"
set xlabel "Number of Haplotypes"
set yrange [0:*]
xmax="`head -1 haplotypes.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'haplotypes.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of # Haplotypes per Contig."
set ylabel "# Haplotypes"
set xlabel "Contigs"
plot 'haplotypes.dat' pt "*" 
pause -1
EOF

#scatter plot of sites vs haplotypes
paste <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f3) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f2) > haplo_sites.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Sites vs Haplotypes per Contig"
set ylabel "# Sites"
set xlabel "# Haplotypes"
plot 'haplo_sites.dat' pt "*" 
pause -1
EOF

#histogram of proportion of individuals haplotyped
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f6 | sort -rg > propindhaplotyped.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of the proportion of individuals haplotyped per contig"
set ylabel "Number of Contigs"
set xlabel "Proportion Haplotyped"
set yrange [0:*]
set xrange [0:1]
binwidth= 0.05
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'propindhaplotyped.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of the proportion of individuals haplotyped per contig"
set ylabel "Proportion Haplotyped"
set xlabel "Contigs"
plot 'propindhaplotyped.dat' pt "*" 
pause -1
EOF

#histogram of number of paralogous individuals
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f8 | sort -rg > paralogs.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of Number of Paralogous Individuals per Contig"
set ylabel "Number of Contigs"
set xlabel "Number of Paralogous Individuals per Contig (0-20)"
set yrange [0:*]
set xrange [0:20]
binwidth= 1
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'paralogs.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of # Paralogous Individuals per Contig"
set ylabel "# Paralogs"
set xlabel "Contigs"
plot 'paralogs.dat' pt "*"
pause -1
EOF

#histogram of number of individuals with low cov/geno errs
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f9 | sort -rg > lowcovgenoerr.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of Number of Individuals With Low Cov/ Geno Errs per Contig"
set ylabel "Number of Contigs"
set xlabel "Number of Individuals per Contig"
set yrange [0:*]
xmax="`head -1 lowcovgenoerr.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'lowcovgenoerr.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Number of Individuals With Low Cov/ Geno Errs per Contig"
set ylabel "Number of Individuals per Contig"
set xlabel "Contigs"
plot 'lowcovgenoerr.dat' pt "*"
pause -1
EOF

#histogram of number of individuals with missing genotypes
tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f10 | sort -rg > missinggeno.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of Number of Individuals With Missing Genotypes per Contig"
set ylabel "Number of Contigs"
set xlabel "Number of Individuals per Contig"
set yrange [0:*]
xmax="`head -1 missinggeno.dat`"
set xrange [0:xmax]
binwidth= xmax/20
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'missinggeno.dat' using (bin($1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Number of Individuals With Missing Genotypes per Contig"
set ylabel "Number of Individuals per Contig"
set xlabel "Contigs"
plot 'missinggeno.dat' pt "*"
pause -1
EOF

#scatter plots of paralogs vs X
paste haplo_sites.dat <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f6) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f8) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f9) <(tail -n +3 $VCF_OUT.Fltr$FILTER_ID.stats.out | cut -f10) > haplo_sites_prophap_paralogs.dat
gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Scatter plot of Paralogous Individuals vs Haplotypes per Contig"
set ylabel "# Paralogs"
set xlabel "# Haplotypes"
plot 'haplo_sites_prophap_paralogs.dat' using 1:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs Sites per Contig"
set ylabel "# Paralogs"
set xlabel "# Sites"
plot 'haplo_sites_prophap_paralogs.dat' using 2:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs Proportion of Individuals Haplotyped per Contig"
set ylabel "# Paralogs"
set xlabel "Proportion Haplotyped"
plot 'haplo_sites_prophap_paralogs.dat' using 3:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs # Individuals with Low Cov/ Geno Err per Contig"
set ylabel "# Paralogs"
set xlabel "# Individuals"
plot 'haplo_sites_prophap_paralogs.dat' using 5:4 with points pt "*"

set title "Scatter plot of Paralogous Individuals vs # Individuals with Missing Genotypes per Contig"
set ylabel "# Paralogs"
set xlabel "# Individuals"
plot 'haplo_sites_prophap_paralogs.dat' using 6:4 with points pt "*"

set title "Scatter plot of # Individuals with Low Cov/ Geno Err vs # Individuals with Missing Genotypes per Contig"
set ylabel "# Individuals w Low Cov / Geno Err"
set xlabel "# Individuals w Missing Geno"
plot 'haplo_sites_prophap_paralogs.dat' using 6:5 with points pt "*"

set title "Scatter plot of # Individuals with Low Cov/ Geno Err vs Proportion of Individuals Haplotyped per Contig"
set ylabel "# Individuals with Low Cov/ Geno Err"
set xlabel "Proportion Individuals Haplotyped"
plot 'haplo_sites_prophap_paralogs.dat' using 3:5 with points pt "*"

set title "Scatter plot of # Individuals with Missing Genotypes vs Proportion of Individuals Haplotyped per Contig"
set ylabel "# Individuals with Missing Genotypes"
set xlabel "Proportion Individuals Haplotyped"
plot 'haplo_sites_prophap_paralogs.dat' using 3:6 with points pt "*"

pause -1
EOF
		
		if [[ PARALLEL == "FALSE" ]]; then
			echo ""; echo -n file $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf has 
			mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			echo SNPs
			echo ""
			echo -n "	Sites remaining:	" && mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | wc -l
			echo -n "	Contigs remaining:	" && mawk '!/#/' $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf | cut -f1 | uniq | wc -l
		else
			echo -n "	Sites remaining:	" && ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz | wc -l " | awk -F: '{a+=$1} END{print a}' 
			echo -n "	Contigs remaining:	" && ls $DataName$CutoffCode.*.bed | parallel --no-notice -k -j $NumProc "tabix -R {} $VCF_OUT.Fltr$FILTER_ID.Haplotyped.vcf.gz | cut -f1 | uniq | wc -l " | awk -F: '{a+=$1} END{print a}'
		fi

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
		
	elif [[ $FILTER_ID == "22" ]]; then
		echo; echo `date` "---------------------------FILTER22: Test For LD -----------------------------"
		
		# grep '^dDocent' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf | cut -f1-2 | tr "_" "\t" | awk 'NR==1 {a1=$1} {printf "%s %.0f\n", $4, ($3*1000)+$4}' | tr " " "\t" | cut -f2 > contig_pos.txt
		# cat <(mawk '/^#/' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf) <(paste <(mawk '!/#/' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf | cut -f1) contig_pos.txt <(mawk '!/#/' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf | cut -f3- ) )> Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.ld.vcf
		# sed -i 's/dDocent_Contig_[0-9]*/fltrVCF_ld_1/g' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.ld.vcf

		grep '^dDocent' $VCF_FILE | cut -f1-2 | tr "_" "\t" | awk 'NR==1 {a1=$1} {printf "%s %.0f\n", $4, ($3*1000)+$4}' | tr " " "\t" | cut -f2 > contig_pos.txt
		cat <(mawk '/^#/' $VCF_FILE) <(paste <(mawk '!/#/' $VCF_FILE | cut -f1) contig_pos.txt <(mawk '!/#/' $VCF_FILE | cut -f3- ) )> ${VCF_FILE%.*}.ld.vcf
		sed -i 's/dDocent_Contig_[0-9]*/fltrVCF_ld_1/g' ${VCF_FILE%.*}.ld.vcf
		Filter="--geno-r2"
		vcftools --vcf ${VCF_FILE%.*}.ld.vcf $Filter
		
	elif [[ $FILTER_ID == "86" ]]; then
		echo; echo `date` "---------------------------FILTER 86: Remove contigs -----------------------------"
		echo "	Before site filter: $VCF_FILE_2"
		echo "	After site filter: $VCF_FILE"
		if [[ $PARALLEL == "TRUE" ]]; then 
			VCF_FILE=${VCF_FILE%.*}
			VCF_FILE_2=${VCF_FILE_2%.*}
			#echo $VCF_FILE $VCF_FILE_2 $VCF_OUT 
			printf "${VCF_FILE}\n${VCF_FILE_2}\n" | parallel --no-notice "gunzip -c {}.gz > {}"
		fi
		VCF_OUT=${VCF_FILE%.*}.Fltr${FILTER_ID}.rmCONTIGS.${COUNTER}.vcf
		echo "	After contig filter: $VCF_OUT"
		comm -23 <(grep '^dDocent_Contig' ${VCF_FILE_2} | sort -g) <(grep '^dDocent_Contig' ${VCF_FILE} | sort -g) | sort -V | cut -f1 | uniq > rmCONTIGS.txt
		grep -v -f "rmCONTIGS.txt" ${VCF_FILE_2} > ${VCF_OUT}
		if [[ $PARALLEL == "TRUE" ]]; then 
			bgzip -@ $NumProc -c ${VCF_OUT} > ${VCF_OUT}.gz
			tabix -p vcf ${VCF_OUT}.gz
		fi
		echo "	$(wc -l rmCONTIGS.txt | cut -f1) contigs filtered"
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

function fltrVCFstats() {
echo -e "File	NumInd	NumContigs	NumSNPs	NumMissingGeno	NumGenoLess10X	NumGeno10-19X	NumGeno20-49X	NumGeno50-99X	NumGeno100-999X" > $DataName$CutoffCode.fltrStats.dat

if [[ $PARALLEL == "TRUE" ]]; then
	ls -tr $DataName$CutoffCode*vcf.gz | parallel -j $NumProc -k "echo -e -n {}'\t' && \  
		zcat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' &&\   #num individuals
		echo -e -n '\t' && \
		zgrep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' && \   #num contigs
		echo -e -n '\t' && \
		zgrep -c '^dDocent' {} | tr -d '\n' && \   #num snps
		echo -e -n '\t' && \
		zgrep -oh '\./\.:' {} | wc -l | tr -d '\n' && \   #num missing genotypes
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[1-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w less than 10x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:1[0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 10-19x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 20-49x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 50-99x cvg
		echo -e -n '\t' && \
		zgrep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l && \   #num genotypes w 100-999x cvg
		" >> $DataName$CutoffCode.fltrStats.dat
else
	ls -tr $DataName$CutoffCode*vcf | parallel -j $NumProc -k "echo -e -n {}'\t' && \  
		cat {} | tail -n 1 | cut -f 10- | tr '\t' '\n' | wc -l | tr -d '\n' &&\   #num individuals
		echo -e -n '\t' && \
		grep '^dDocent' {} | cut -f1 | uniq | wc -l | tr -d '\n' && \   #num contigs
		echo -e -n '\t' && \
		grep -c '^dDocent' {} | tr -d '\n' && \   #num snps
		echo -e -n '\t' && \
		grep -oh '\./\.:' {} | wc -l | tr -d '\n' && \   #num missing genotypes
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[1-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w less than 10x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:1[0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 10-19x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[2-4][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 20-49x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[5-9][0-9]:' {} | wc -l | tr -d '\n' && \   #num genotypes w 50-99x cvg
		echo -e -n '\t' && \
		grep -oh '[01]/[01]:[1-9][0-9][0-9]:' {} | wc -l && \   #num genotypes w 100-999x cvg
		" >> $DataName$CutoffCode.fltrStats.dat
fi
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

		[summary stats]
                -S				calculate the number of individuals, contigs, SNPs, and bin the
								 genotypes by depth of coverage for each filter step and compile
								 into one output file *.fltrVCFstats.dat
		

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
echo "	R"
echo "		tidyverse"
echo "		gridExtra"
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
while getopts ":f:c:b:d:v:g:p:s:w:r:o:t:PSh" opt; do
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
	S)
        echo "	Forcing all filters to run in parallel:	$OPTARG"
		STATS=TRUE
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

if [[ $STATS == "TRUE" ]]; then
	echo "	-S set. Summary filter data will be saved to a *fltrVCFstats.dat file."
else
	STATS=FALSE
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

if [[ $STATS == "TRUE" ]]; then
	echo ""; echo `date` " --------------------------- Compile summary stats! ---------------------------"; echo "" 
	fltrVCFstats
	echo ""; echo `date` " --------------------------- Summary stats complete! ---------------------------"; echo "" 
fi