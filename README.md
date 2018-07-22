# IMPACT (Inference and Modeling of Phenotype-related ACtive Transcription) 

A data aggegation strategy to model specific transcriptional processes based on transcription factor and cell-state specificity. 

## Biorxiv preprint
```
https://www.biorxiv.org/content/early/2018/07/10/366864
```

## Using IMPACT 

### Directory Contents: 


Data: 

Motif/ Transcription factor binding motif PWMs (position weight matrices) and genome-wide positions of motifs 

Features/ Epigenomic and sequence features used in the IMPACT framework to predict cell-state-specific regulatory elements. 


Implementing IMPACT: 

Training/ use Pipeline.sh to create positive (regulatory) and negative (non-regulatory) training sets. 

GenomeWide/IMPACT_modelfit/ use ENet_fit7_TFArgument.R to run model on a per-TF basis, to get betas on epigenomic and sequence features 

GenomeWide/GenomeTracks/Predictvalues_TFarg.R predicts cell-state-specific regulatory elements genomewide based on model fits 

GenomeWide/GenomeTracks/MinMaxValues_TFarg.R computes minimum and maximum values per chromosome

GenomeWide/GenomeTracks/Scale_TFarg.R scales logistic regression output to be between 0 and 1. 

GenomeWide/ldsc/ use CreateAnnotations_TFarg.R to create IMPACT annotations formatted to be run through S-LDSC to compute heritability estimates

sLDSC/ executables to run S-LDSC in European and East Asian populations. 


Assessing IMPACT performance: 

Performance/ compute AUC and sensitivity of IMPACT to predict TF binding of a motif.

KnockOuts/ compute AUC and sensitivity after parameter variation and categorical feature knock outs.

CandidateLoci/ find loci with predicted IMPACT regulatory elements that do not have evidence of TF binding.

TargetvNonTargets/ compare IMPACT predictions at TF target genes, defined by ChIP-seq, to TF non-target genes.

