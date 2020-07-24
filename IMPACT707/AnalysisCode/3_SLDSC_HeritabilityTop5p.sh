module load gcc/6.2.0
module load R/3.3.3

annotnum=655 #for example 
trait=PASS_Rheumatoid_Arthritis #for example
shortname=RA #for example

annot=IMPACT_ENCODE_${annotnum}
annot_dir=/path/to/annots
quantiles=20
pop=EUR
cmd="Rscript /path/to/ldsc/quantile_h2g.r custBLcts.$annotnum.q$quantiles.${pop}.newM /path/to/Results/${pop}/Kawakami/$annot/$trait/$annot.ldsc.custBLcts custBLcts.$annot.$shortname.q$quantiles.${pop}.new.txt" 
bsub -q big -R 'rusage[mem=15000]' -W 10:00 -o $annot.log "$cmd"

