Scan genome for motifs in this directory.

Replace Motif_TF_cellstate.txt with your TF and cellstate. 

Run in bash. 

perl /data/srlab/amariuta/homer/bin/scanMotifGenomeWide.pl Motif_TF_cellstate.txt hg19 > scanMotifsgenomewide.TF_cellstate.txt
sort -t $'\t' -k6,6rn scanMotifsgenomewide.TF_cellstate.txt > scanMotifsgenomewide.TF_cellstate.sort.txt
head -n 15000 scanMotifsgenomewide.TF_cellstate.sort.txt > scanMotifsgenomewide.TF_cellstate.15000.sort.txt
