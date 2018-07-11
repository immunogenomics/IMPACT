create_annot=1
trait_file=/home/tja2/TFBStrack/Tbet/RA.trait.path.txt

annot_main_path=/home/tja2/TFBStrack/Amariutaetal/
annots=(Tbet Gata3 Stat3 Foxp3)
for tf in ${annots[@]}
do

	if [ $create_annot -eq 1 ]
	    then
	    for chrom in {1..22}
	    do

		echo "Working on chr$chrom"
		sbatch -p short -t 0-10:00 --mem=15000 --wrap="python /n/groups/price/bryce/ldsc/ldsc.py --bfile /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --out $annot_main_path/$tf/${tf}_GenomewideTrack_IMPACT.$chrom --yes-really --annot $annot_main_path/$tf/${tf}_GenomewideTrack_IMPACT.$chrom.annot.gz --print-snps /home/fh80/Code/RunLDSC/list.txt"
	    done
	fi


	if [ $create_annot -eq 2 ]
	then
	    for trait in `cat $trait_file`
	    do
		t=`basename $trait | sed 's/.sumstats//g'`
		mkdir $annot_main_path/$tf/Results
		mkdir $annot_main_path/$tf/Results/${tf}_GenomewideTrack_IMPACT
		mkdir $annot_main_path/$tf/Results/${tf}_GenomewideTrack_IMPACT/$t
		echo "Working on trait $t"
		sbatch -p short -t 0-10:00 --mem=15000 --wrap="python /n/groups/price/steven/soft/ldsc/ldsc.py --h2 $trait --ref-ld-chr /home/tja2/TFBStrack/Amariutaetal/BaselineLD_customized/customized/customized_baselineLD_cts.,$annot_main_path/$tf/${tf}_GenomewideTrack_IMPACT. --w-ld-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC. --overlap-annot --frqfile-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC. --out $annot_main_path/$tf/Results/${tf}_GenomewideTrack_IMPACT/$t/${tf}_GenomewideTrack_IMPACT.${t}.ldsc.custBLcts --print-coefficients --print-delete-vals"
	    done
	fi
done


