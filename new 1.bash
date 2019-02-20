		grep '^dDocent' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf | cut -f1-2 | tr "_" "\t" | awk 'NR==1 {a1=$1} {printf "%s %.0f\n", $4, ($3*1000)+$4}' | tr " " "\t" | cut -f2 > contig_pos.txt
		cat <(mawk '/^#/' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf) <(paste <(mawk '!/#/' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf | cut -f1) contig_pos.txt <(mawk '!/#/' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.vcf | cut -f3- ) )> Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.ld.vcf
		sed -i 's/dDocent_Contig_[0-9]*/fltrVCF_ld_1/g' Hlobatus.D2.5.5.Fltr21.1.MostInformativeSNP.ld.vcf
