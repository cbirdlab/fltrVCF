
#calculate mean of site mean depths by contig
awk '{x[$1] += $3; N[$1]++} END{for (i in x) print i, N[i], x[i]/N[i]}' PIRE_SiganusSpinus.D.15.15.Fltr04.2.ldepth.mean.150 | tr -s " " "\t" | sort -k3 -n | less -S

awk '{x[$1] += $3; N[$1]++} END{for (i in x) print i, N[i], x[i]/N[i]}' PIRE_SiganusSpinus.D.15.15.Fltr04.2.ldepth.mean.150 | tr -s " " "\t" | sort -k3 -n | cut -f3 | awk '{all[NR] = $0 } END{print all[int(NR*0.99 - 0.5)]}'

