
usage 1: 
bsub < ENet_fit7_TFArgument.lsf 
set TF=yourTF in ENet_fit7_TFArgument.lsf

usage 2:
interactive mode 
Rscript ENet_fit7.R [arg1 1st TF to model] [arg2 2nd TF to model] ... 

example: 

Rscript ENet_fit7.R Tbet Gata3 Stat3 Foxp3


