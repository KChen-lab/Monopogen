

library(e1071)
library(reshape2)

source("/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/src/svm_phasing_sub.R")

args = commandArgs(trailingOnly=TRUE)

chr <- args[1]
loaddata <- args[2]
dis_limit_physical <- as.numeric(args[3])
dis_limit_genetic <- as.numeric(args[4])





print(chr)

if(loaddata==1){
	print("Loading data now...")
	mat <- read.table(paste0(file="../germline/chr",chr,".gl.filter.hc.cell.mat.gz"))
	#normal <- read.table(file=paste0("normal.",chr,".snv"))
	#tumor <- read.table(file=paste0("tumor.",chr,".snv"))
	#somatic <- read.table(file="somatic.snv")
	feature_info <- read.csv(paste0(file="../debug/chr",chr,".gl.feature.info"),header=F)


	mat$V2 <- paste0(mat$V1,":",mat$V2)
	meta <- mat[,c("V1","V2","V3","V4")]
	meta <- as.data.frame(meta)
	meta$in_normal <- 0
	meta$in_tumor <- 0
	meta$somatic <- 0
	#meta$in_normal[meta$V2%in%normal$V1] <- 1
	#meta$in_tumor[meta$V2%in%tumor$V1] <- 1
	#meta$somatic[meta$V2%in%somatic$V1] <- 1
	meta$baf <- frq(mat)

	#meta$baf <-0.5
	# SVM module on remove low quality mutations


	svm_info <- feature_info[feature_info$V1%in%meta$V2,]
	svm_dt <- getFeatureInfo(svm_info=svm_info)

	overlap <- intersect(meta$V2, svm_dt$id)

	meta_qc <- meta[meta$V2%in%overlap,]
	mat_qc <- mat[meta$V2%in%overlap,]
	svm_dt_qc <- svm_dt[svm_dt$id%in%overlap,]
	svm_dt_qc$BAF <- meta_qc$baf
	save.image(file=paste0(chr,".input.image.RDS"))
}

if(loaddata==2){
	load(file=paste0(chr,".input.image.RDS"))
}


source("/rsrch3/scratch/bcb/jdou1/scAncestry/Monopogen/src/svm_phasing_sub.R")
print(args)


print("start to run SVM and haplotype filtering")
print("new version ")

# get putative mutation block 
#meta_qc$flag <- paste0(meta_qc$in_normal, ":", meta_qc$in_tumor, ":", meta_qc$somatic)





mutation_block <- SNV_block(summary=meta_qc)

svm_in <- SVM_prepare(mutation_block)

svm_out <- SVM_train(label =svm_in, data=svm_dt_qc, downsampling=2)

p_lower <- sort(svm_out$test$POS)[floor(nrow(svm_out$test)*0.1)]
p_upper <- sort(svm_out$test$POS)[floor(nrow(svm_out$test)*0.95)]
filter1 <- (svm_out$test$POS>p_lower)
filter2 <- (svm_out$test$baf>0.1 )
svm_pass <- svm_out$test[filter1 & filter2,]


svm_phasing <- Phasing(x=mat_qc, snv_meta=svm_pass, dis_limit=dis_limit_genetic, readlength=dis_limit_physical)

svm_out$phasing <- svm_phasing
svm_out$feature <- svm_dt_qc

saveRDS(svm_out,file=paste0("svm",chr,".out.RDS"))



x=mat_qc
snv_meta=svm_pass
dis_limit=dis_limit_genetic
readlength=dis_limit_physical