library(polyester)
library(Biostrings)

fasta<-readDNAStringSet("./chr22.fa")

set.seed(10)

#Create a vector for the FC of each gene between groups, this parameter is mandatory
case<-sample(c(4/1, 2/1, 3/2,1/1,2/3,1/2,1/4),
	size=(length(fasta)),
	replace=TRUE,
	prob=c(0.03,0.07,0.09,0.7,0.06,0.03,0.02))

#Create a FC matrix for the control group
fold_chages<-as.matrix(data.frame(case, "control"=1))
colnames(fold_chages)<-NULL

readspert<-round(20 * width(fasta)/100)

simulate_experiment("./chr22.fa", reads_per_transcript=readspert,
	error_model="illumina4",
	num_reps=c(12,12), fold_changes=fold_chages, seed=10,
	out_dir="sim_baixa_var",
	)
	
simulate_experiment("./chr22.fa", reads_per_transcript=readspert,
	error_model="illumina4",
	num_reps=c(12,12), fold_changes=fold_chages, seed=10,
	out_dir="sim_alta_var",
	size=1 #Parameter for high variance model
	)
