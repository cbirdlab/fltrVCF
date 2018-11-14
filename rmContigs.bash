#!/bin/bash

#This is a script that identifies which contigs had SNPs that were filtered, then filters those contigs

#grep -Fxv -f dDocent_Contig_4021_removed_pos_122.vcf dDocent_Contig_4021.vcf

#comm -23 dDocent_Contig_4021.vcf dDocent_Contig_4021_removed_pos_122.vcf


#in the next line, the first file is the orig and the second is filtered
RMcontigs=($(comm -23 dDocent_Contig_4021.vcf dDocent_Contig_4021_removed_pos_122.vcf | cut -f1))
echo ${RMcontigs[@]} | parallel --no-notice -k -j ${NumProc} "grep -v ${RMcontigs} dDocent_Contig_4021.vcf" 


