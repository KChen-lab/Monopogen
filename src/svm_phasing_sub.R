
frq <- function(x=mat) {
	res<-matrix(0,nrow(mat),1)
	for(i in seq(1,nrow(mat),1)){
		geno <- mat[i, seq(5,ncol(mat),1)]
		geno <- gsub("/","|", geno)
		geno <- geno[!geno=="0|0"]
		s <- strsplit(geno, "|")
		result <- sapply(s, "[[", 1)
		if(length(result)>10){

			myfrq <- length( which(result=="0"))/length(result)
		}else{
			myfrq <- 0
		}
		if(myfrq>0.5){myfrq <- 1- myfrq}
		res[i,1] <- myfrq
	}
	return(res)
}


getFeatureInfo <-function(svm_info=NULL){
	n_min <-(-1000)
	dt <- data.frame("id"=svm_info$V1, "DP"=n_min, "QS"=n_min, 
		"VDB"=n_min,"SGB"=n_min,"RPB"=n_min, "MQB"=n_min,"MQSB"=n_min, "BQB"=n_min)

	svm_info[is.na(svm_info)] <-"-1000"
	for(i in seq(1,nrow(svm_info),1)){
		for(j in seq(1,ncol(svm_info),1)){
			if(svm_info[i,j]=="DP"){dt$DP[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="QS"){dt$QS[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="VDB"){dt$VDB[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="SGB"){dt$SGB[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="RPB"){dt$RPB[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="MQB"){dt$MQB[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="MQSB"){dt$MQSB[i] <- svm_info[i,j+1]}
			if(svm_info[i,j]=="BQB"){dt$BQB[i] <- svm_info[i,j+1]}
		}
	}
	return(dt)
}

SNV_block <-function(summary=NULL){
	block <- 1
	last <- ".|."
	summary$region <- 0
	for(i in seq(1,nrow(summary),1)){
		if(summary[i,4]==".|."){
			if(last==".|." | last=="0|0"){;}
			else{	
				block <- block + 1
			}
			summary$region[i] <- block
		}
		last <- summary[i,4]
	}
	return(summary)
}


SVM_prepare <-function(x=NULL){
	svm<-list()
	a<- table(x$region)
	neg <- x[(x$region%in% names(a)[a>=4] & x$V4==".|.") | (x$V4=="0|0" & x$baf>0),]
	#neg <- x[(x$region%in% names(a)[a>=3] & x$V4==".|.") | x$V4=="0|0",]
	#neg <- x[x$V4=="0|0" | ,]
	pos <- x[!x$V4==".|." & !x$V4=="0|0" &x$baf>0.4,]
	svm$neg <- neg 
	svm$pos <- pos
	svm$test <- x[x$region%in% names(a)[a<4] & x$V4==".|.",]
	return(svm)
}

SVM_train <- function(label=NULL, data=NULL, downsampling=20){
	

	print("Run SVM trainning")

	feature <-c("QS", "VDB", "SGB", "RPB", "MQB", "MQSB", "BQB", "BAF")
	data$BQB <- as.numeric(data$BQB)

	data$VDB[data$VDB==(-1000)] <-0
	data$MQSB[data$MQSB==(-1000)] <-1
	data$QS <- abs(data$QS-0.5)
	x_pos <- data[data$id%in%label$pos$V2,feature]
	#x_pos <- x_pos[seq(1,nrow(x_pos),downsampling),]
	x_neg <- data[data$id%in%label$neg$V2,feature]

	y_pos <- rep("POS", nrow(x_pos))
	y_neg <- rep("NEG", nrow(x_neg))

	if(nrow(x_pos) > nrow(x_neg)){
		# randomly select matched variants between pos and neg
		sel <- sample(nrow(x_pos))[seq(1,nrow(x_neg),1)]
		x_pos <- x_pos[sel,]
		y_pos <- y_pos[sel]
	}
	if(nrow(x_pos) < nrow(x_neg)){
		# randomly select matched variants between pos and neg
		sel <- sample(nrow(x_neg))[seq(1,nrow(x_pos),1)]
		x_neg <- x_neg[sel,]
		y_neg <- y_neg[sel]
	}
    
	print(paste0("pos->:", length(y_pos)))
	print(paste0("neg->:", length(y_neg)))

	model <- svm(as.matrix(rbind(x_pos, x_neg)), as.factor(c(y_pos,y_neg)), probability=TRUE)
	test <- data[data$id%in%label$test$V2,feature]
	pred <- predict(model, as.matrix(test), probability=TRUE)
	prob <- attr(pred, "probabilities")
	prob <-as.data.frame(prob)
	label$test$POS <- prob$POS
	label$test$NEG <- prob$NEG
	return(label)
}

calcP <- function(geno_vec, dis_vec, dis_limit=20000){

		dis_bin <-  c(100,  1000,  5000,   10000, 20000,  500000, 1000000000000000000)
		prob_bin <- c(1,    0.93,  0.78,   0.71,  0.69,   0.64,   0)

		prob_bin[dis_bin>dis_limit] <- 0

		pos_match <- which(geno_vec%in%c("000","111"))
		dis_match <- dis_vec[pos_match]

		match_p <- 0
		for(pos in dis_match){
				if(pos<=dis_bin[1]){p <-1}
				else{
					p <- prob_bin[max(which(pos>=dis_bin))+1]
				}
				match_p <- match_p + p/sqrt(pos) 
		}

		pos_mismatch <- which(geno_vec%in%c("010","101"))
		dis_mismatch <- dis_vec[pos_mismatch]


		mismatch_p <- 0
		for(pos in dis_mismatch){
				if(pos<=dis_bin[1]){p <-1}
				else{
					p <- prob_bin[max(which(pos>=dis_bin))+1]
				}
				mismatch_p <-mismatch_p + p/sqrt(pos) 
		}

		if((match_p + mismatch_p)==0){p <- (-1)}
		else{
			p <- match_p/(match_p+mismatch_p)
		}
		
		return(p)
}

Phasing <- function(x = NULL, snv_meta = NULL, dis_limit = 20000, readlength=NULL){

	start <- 1
    if(start){
		y <- x
	    dis <- colsplit(y$V1,":", c("id1","id2"))[,2]
		dis_rd_study <-c()
		trio_hap_rd_study <-c()
		hap_rd2 <-c()
		dis_rd2 <-c()
		cnt <- 0
	
		check_study <- snv_meta
		check_study$prob <- (-1000)
		check_study$tol <-  0
		check_study$dis <- (-1000)

		check_study$prob_phy <- (-1000)
		check_study$tol_phy <- 0
		res<-c()

		pos_lst_study <-list()
		for(i in seq(5,ncol(y),1)){
		 	vec <- y[,i]
		 	pos_lst_study[[i]] <- which(!vec=="0|0" & !y$V4==".|." & !vec=="0/0")
		}
	}

	for(snv in snv_meta$V2){
	 	  cnt=cnt+1
	 	  #print(paste0("scanning ",cnt," now..."))
	 	  snv_index <- which(y$V2==snv)
	 	  snv_index_noMisCol <- which(!y[snv_index,]=="0|0" & !y[snv_index,]=="0/0",)

	 	  # remove the first column since they are not genotypes 
	 	  snv_index_noMisCol <- snv_index_noMisCol[snv_index_noMisCol>4]
	 	  geno_vec <-c()
	 	  dis_vec <-c()
	 	  dis_phy1 <- c()
		  dis_phy2 <- c()
	 	  # identify haplotypes for each element
	 	  index <-0
	 	  for(pos in snv_index_noMisCol){
	 	  	 # element index (snv_index, pos)
	 	  	 pos_noMisCol <- pos_lst_study[[pos]]
	 	  	 lower_vec <- which(pos_noMisCol<snv_index)
	 	  	 upper_vec <- which(pos_noMisCol>snv_index) 	  	 
	 	  	 if(length(lower_vec)>=1 & length(upper_vec)>=1){
	 	  	 		lower <- pos_noMisCol[max(lower_vec)]
	 	  	 		upper <- pos_noMisCol[min(upper_vec)]
	 	  	 		tp <- y[snv_index, pos]
	 	  	 		tp <- gsub("/", "|", tp) 
	 	  	 		geno <- c(y[lower,pos], tp, y[upper, pos])

	 	  	 		part <-colsplit(geno,"\\|",c("i1","i2"))[,2]
					part[part>1] <- 1
					flag <- paste0(part[1],part[2],part[3])
					hap_rd2 <-c(hap_rd2, abs(part[2]-part[1]), abs(part[3]-part[2]))
					dis_rd2 <-c(dis_rd2, dis[snv_index]-dis[lower], dis[upper]-dis[snv_index])

	 	  	 		dis_vec <- c(dis_vec, dis[upper]-dis[lower])
	 	  	 		dis_phy1 <- c(dis_phy1,  dis[snv_index]-dis[lower])
	 	  	 		dis_phy2 <- c(dis_phy2, dis[upper]-dis[snv_index])
	 	  	 		geno_vec <-c(geno_vec,flag)
	 	  	 }
	 	  }

	 	  sum <- 0
	 	  sum_equal <- 0
	 	  if(sum(dis_phy1<readlength)){
	 	  	 d_geno <- substr(geno_vec[which(dis_phy1<readlength)],1,2)
	 	  	 sum <- length(d_geno)
	 	  	 sum_equal <- sum(d_geno=="00" | d_geno=="11")
	 	  }

	 	  if(sum(dis_phy2<readlength)){
	 	  	 d_geno <- substr(geno_vec[which(dis_phy2<readlength)],2,3)
	 	  	 sum <- length(d_geno) + sum
	 	  	 sum_equal <- sum(d_geno=="00" | d_geno=="11") + sum_equal
	 	  }


	 	  if(sum>0){
	 	  	  		p <- sum_equal/sum
	 	  	  		check_study$prob_phy[cnt] <- p
	 	  	  		check_study$tol_phy[cnt] <- sum
	 	  	  		print(paste0("Scan ", cnt,":", snv, " ", sum, " ", p))
	 	  }

	 	  if(length(geno_vec)>0 ){
	 	  	  dis_rd_study <-c(dis_rd_study, dis_vec)
			  trio_hap_rd_study <-c(trio_hap_rd_study, geno_vec)
	 	  	  #print(table(trio_hap_rd))
	  	  		check_study$prob[cnt] <- calcP(geno_vec,dis_vec,dis_limit)
	  	  		check_study$tol[cnt] <- sum(dis_vec<dis_limit)
	  	  		#check_study$dis[cnt] <- max(dis_vec[dis_vec<dis_limit])
	 	  }
	}
	res$table <- check_study
	res$dis_genetic <- dis_rd_study 
	res$hap_genetic <- trio_hap_rd_study
	res$dis_physical <- dis_rd2 
	res$dis_genetic <- hap_rd2

	return(res)
}




