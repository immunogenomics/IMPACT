
#Get S-LDSC results together
pop=EAS #or EUR

for i in {1..728}
do
        cd /path/to/Results/${pop}/Kawakami/IMPACT_ENCODE_${i}/
        for j in *
        do
        cd /path/to/Results/${pop}/Kawakami/IMPACT_ENCODE_${i}/${j} 
        checkfile=IMPACT_ENCODE_${i}.ldsc.custBLcts.results
        cp ${checkfile} /path/to/datadump/IMPACT_ENCODE_${i}.${j}.${pop}.ldsc.custBLcts.results
done

cd /path/to/datadump
for i in *
do
tail -n 1 $i >> IMPACTxTRAIT_KawakamiIMPACT_${pop}.txt
echo $i >> IMPACTxTRAIT_KawakamiIMPACT_${pop}_filekey.txt
done
gzip /path/to/datadump/IMPACTxTRAIT_KawakamiIMPACT_${pop}.txt
gzip /path/to/datadump/IMPACTxTRAIT_KawakamiIMPACT_${pop}_filekey.txt

#Get trait heritability estimates together 
pop=EUR
for i in {1..728}
do
        cd /path/to/Results/${pop}/Kawakami/IMPACT_ENCODE_${i}/
        for j in *
        do
        cd /path/to/Results/${pop}/Kawakami/IMPACT_ENCODE_${i}/${j}
	checkfile=IMPACT_ENCODE_${i}.ldsc.custBLcts.log
        grep Observed ${checkfile} >> /path/to/datadump/IMPACTxTRAIT_KawakamiIMPACT_${pop}_h2.txt
        echo $i,$j >> /path/to/datadump/IMPACTxTRAIT_KawakamiIMPACT_${pop}_h2_key.txt
done
done
gzip /path/to/datadump/IMPACTxTRAIT_KawakamiIMPACT_${pop}_h2.txt
gzip /path/to/datadump/IMPACTxTRAIT_KawakamiIMPACT_${pop}_h2_key.txt

