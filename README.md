# IMPACT (Inference and Modeling of Phenotype-related ACtive Transcription) 

A data aggregation strategy to model specific transcriptional processes based on transcription factor and cell-state specificity. 

## Publication in American Journal of Human Genetics: IMPACT: Genomic Annotation of Cell-State-Specific Regulatory Elements Inferred from the Epigenome of Bound Transcription Factors
* [AJHG](https://www.cell.com/ajhg/supplemental/S0002-9297%2819%2930108-9) 

## Publication in Nature Genetics: Improving the trans-ancestry portability of polygenic risk scores by prioritizing variants in predicted cell-type-specific regulatory elements
* [Nature Genetics](https://www.nature.com/articles/s41588-020-00740-8)

## Using IMPACT 


### Dependencies 

1. A motif finding algorithm, such as HOMER as used in our implementation of IMPACT. 
2. R version 3.2.2 or higher
3. R package GenomicRanges

### Data 

Transcription factor binding motif PWMs (position weight matrices) and genome-wide positions of motifs 
```
Motifs/ 
```

Epigenomic and sequence features used in the IMPACT framework to predict cell-state-specific regulatory elements. 
```
Features/ 
```

### Implementing IMPACT 
1. Create positive (regulatory) and negative (non-regulatory) training sets:
```
Training/ Pipeline.sh "TF" "cell type" "number input (ChIP-seq) files"
```

2. Run model on a per-TF basis, to get betas on epigenomic and sequence features:
```
GenomeWide/IMPACT_modelfit/ENet_fit7_TFArgument.R "TF"
```

3. Predict cell-state-specific regulatory elements genomewide based on model fits: 
```
GenomeWide/GenomeTracks/Predictvalues_TFarg.R "TF"
```

4. Then compute minimum and maximum values per chromosome:
```
GenomeWide/GenomeTracks/MinMaxValues_TFarg.R "TF"
```

5. Finally, scale logistic regression output to be between 0 and 1: 
```
GenomeWide/GenomeTracks/Scale_TFarg.R "TF"
```

6. Create IMPACT annotations formatted to be run through S-LDSC to compute heritability estimates:
```
GenomeWide/ldsc/CreateAnnotations_TFarg.R "TF"
```

7. Run S-LDSC in European and East Asian populations: 
```
sLDSC/ 
```

### Assessing IMPACT performance: 

Compute AUC and sensitivity of IMPACT to predict TF binding of a motif:
```
Performance/ 
```
Compute AUC and sensitivity after parameter variation and categorical feature knock outs:
```
KnockOuts/ 
```

Find loci with predicted IMPACT regulatory elements that do not have evidence of TF binding:
```
CandidateLoci/ 
```

Compare IMPACT predictions at TF target genes, defined by ChIP-seq, to TF non-target genes.
```
TargetvNonTargets/ 
```


