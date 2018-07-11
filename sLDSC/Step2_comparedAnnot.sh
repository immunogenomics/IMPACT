trait_file=/home/tja2/TFBStrack/Tbet/RA.trait.path.txt
annot=Tbet

            for trait in `cat $trait_file`
            do
                t=`basename $trait | sed 's/.sumstats//g'`
                mkdir /home/tja2/TFBStrack/Amariutaetal/ComparisonToStandardAnnotations/comparedAnnot
                mkdir /home/tja2/TFBStrack/Amariutaetal/ComparisonToStandardAnnotations/comparedAnnot/$annot
                mkdir /home/tja2/TFBStrack/Amariutaetal/ComparisonToStandardAnnotations/comparedAnnot/$annot/$t
                echo "Working on trait $t"
                sbatch -p short -t 0-10:00 --mem=15000 --wrap="python /n/groups/price/steven/soft/ldsc/ldsc.py --h2 $trait --ref-ld-chr /home/tja2/TFBStrack/Amariutaetal/BaselineLD_customized/customized/customized_baselineLD_cts.,/home/tja2/TFBStrack/Amariutaetal/$annot/${annot}_GenomewideTrack_IMPACT.,/home/tja2/TFBStrack/Amariutaetal/comparedAnnot/comparedAnnot. --w-ld-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC. --overlap-annot --frqfile-chr /n/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC. --out /home/tja2/TFBStrack/Amariutaetal/ComparisonToStandardAnnotations/comparedAnnot/comparedAnnot/$t/comparedAnnot.${t}.ldsc.custBLcts --print-coefficients --print-delete-vals"
            done

