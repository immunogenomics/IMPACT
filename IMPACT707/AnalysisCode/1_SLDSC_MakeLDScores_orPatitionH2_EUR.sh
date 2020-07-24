#README: for all /path/to/? variants, please specify your own directories accordingly 
trait_file=/path/to/sumstats/sumstats.txt #each line is a path to sumstats, e.g. 69 lines for 69 EUR sumstats
create_annot=1 #if 1, make LD scores, if 2 partition h2 using sumstats

for file in {1..728}
do
	annot=IMPACT_ENCODE_${file}
	annot_dir=/path/to/annots #specify this using your directory
	if [ $create_annot -eq 1 ]
	    then
	    for chrom in {1..22}
	    do

		echo "Working on chr$chrom"
		cmd="python /path/to/ldsc/ldsc.py --bfile /path/to/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chrom --l2 --ld-wind-cm 1 --out $annot_dir/$annot.$chrom --yes-really --annot $annot_dir/$annot.$chrom.annot.gz --print-snps /path/to/ldsc/SNP_rsID_list_EUR.txt"
                bsub -q big -R 'rusage[mem=15000]' -W 10:00 -o $chrom.$annot.log "$cmd"
	    done
	fi
    
        if [ $create_annot -eq 2 ]
        then
            for trait in `cat $trait_file`
            do
                t=`basename $trait | sed 's/.sumstats//g'`
		mkdir Results/EUR/Kawakami
                mkdir Results/EUR/Kawakami/$annot
                mkdir Results/EUR/Kawakami/$annot/$t
                echo "Working on trait $t"
		cmd="python /path/to/ldsc/ldsc.py --h2 $trait --ref-ld-chr /path/to/ldsc/customized_baselineLD_cts.,$annot_dir/$annot. --w-ld-chr /path/to/ldsc/weights.hm3_noMHC. --overlap-annot --frqfile-chr /path/to/ldsc/1000G_Phase3_frq/1000G.EUR.QC. --out Results/EUR/Kawakami/$annot/$t/${annot}.ldsc.custBLcts --print-coefficients --print-delete-vals"
                bsub -q big -R 'rusage[mem=15000]' -W 10:00 -o $chrom.$annot.log "$cmd"
            done
        fi


done
