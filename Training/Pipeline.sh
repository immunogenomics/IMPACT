echo "First arg: $1" #TF
echo "Second arg: $2" #cell-type
echo "Third arg: $3" #number of input files

module load R/3.2.2
scanGenome=2 #1- need to scan genome for motifs, 2- files already made
createannotations=2 #1- need to create annotations, 2-don't need to
makemotifs=1 #1- need to make .motif files, 2- don't need to
#in current IMPACT, we only search for 1 motif per annotation. script may handle > 1 motif (for multiple TF IMPACT annotations). 

Rscript Step1_InputChIPandFeatures.R $1 $2 $3

if [ $makemotifs -eq 1 ]
	then
	num_motifs=$(ls -d *TF.* | wc -l)
	for i in $(seq 1 $num_motifs);
	do
	dos2unix TF.$i.PWM.txt
	cp TF.$i.PWM.txt TF.$i.PWM.motif
	rm TF.$i.PWM.txt
	done
	dos2unix AllTFs_PWM.txt
	cp AllTFs_PWM.txt AllTFs_PWM.motif
	rm AllTFs_PWM.txt
fi

#homer find instances of motifs
perl /data/srlab/amariuta/homer/bin/findMotifsGenome.pl train_df.txt hg19 homer_output/ -size given -find AllTFs_PWM.motif > TFs.findinstances.txt

#build positive set
#here, change training sample size
Rscript Step3_ChIPwithMotifs_PositiveSet_1motif_TTtogether_update.R

#negative set
num_lines=15000

if [ $scanGenome -eq 1 ]
	then
	for i in $(seq 1 $num_motifs); 
	do
	perl /data/srlab/amariuta/homer/bin/scanMotifGenomeWide.pl TF.$i.PWM.motif hg19 > scanMotifsgenomewide.$i.txt
	sort -t $'\t' -k6,6rn scanMotifsgenomewide.$i.txt > scanMotifsgenomewide.$i.sort.txt	
	head -n $num_lines scanMotifsgenomewide.$i.sort.txt > scanMotifsgenomewide.$i.$num_lines.sort.txt
	done
fi

if [ $scanGenome -eq 2 ]
        then
	cp scanMotifsgenomewide.${1}_${2}.15000.sort.txt scanMotifsgenomewide.1.${num_lines}.sort.txt
fi

Rscript Step4_NegativeSet_1motif_TTtogether_update.R

cp train_test_positive_bed.txt train_test_positive_bed_${1}only_center.txt 
cp train_test_negative_bed.txt train_test_negative_bed_${1}only_center.txt

Rscript PhastconsConservation_center.R ${1}














