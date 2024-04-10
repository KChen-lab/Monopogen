library(reshape2)
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

		if(is.na(match_p) | is.na(mismatch_p)){p<-(-1)}
		else{
			if((match_p + mismatch_p) ==0){
				p <- (-1)
			}else{
				p <- match_p/(match_p+mismatch_p)
			}
		}
		
		return(p)
}
	


  

args = commandArgs(trailingOnly=TRUE)

chr <- args[1]
out <- args[2]
len <- as.numeric(args[4])
bin <- as.numeric(args[5])

load(paste0(args[3],chr,".input.image.RDS"))



germline <- mat_qc[!mat_qc$V4=="0|0" & !mat_qc$V4==".|.",]

print(dim(germline))

y <- germline
dis <- colsplit(y$V1,":", c("id1","id2"))[,2]
check <- germline[,c(1,2,3,4)]
check$prob <- (-1000)
check$tol <- 0
check$dis <- (-1000)

check$prob_phy <- (-1000)
check$tol_phy <- 0

pos_lst <-list()
sta <- 0
for(i in seq(5,ncol(y),1)){
 	vec <- y[,i]
 	pos_lst[[i]] <- which(!vec=="0|0" & !y$V4==".|.")
 	sta <- sta + length(which(!vec=="0|0" & !y$V4==".|."))
}


dis_rd <-c()
trio_hap_rd <-c()
cnt <- 0
dis_limit <- 20000

lst <- germline$V2
dis_rd2 <-c()
hap_rd2 <-c()

for(snv in lst[seq(1,length(lst),bin)]){
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
 	  	 pos_noMisCol <- pos_lst[[pos]]
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

 	  if(!is.null(dis_phy1) & sum(dis_phy1<len)){
 	  	 d_geno <- substr(geno_vec[which(dis_phy1<len)],1,2)
 	  	 sum <- length(d_geno)
 	  	 sum_equal <- sum(d_geno=="00" | d_geno=="11")
 	  }

 	  if(!is.null(dis_phy2) & sum(dis_phy2<len)){
 	  	 d_geno <- substr(geno_vec[which(dis_phy2<len)],2,3)
 	  	 sum <- length(d_geno) + sum
 	  	 sum_equal <- sum(d_geno=="00" | d_geno=="11") + sum_equal
 	  }


 	  if(sum>0){
 	  	  		p <- sum_equal/sum
 	  	  		check$prob_phy[cnt] <- p
 	  	  		check$tol_phy[cnt] <- sum
 	  	  		print(paste0("Scan ", cnt,":", snv, " ", sum, " ", p, " ",len))
 	  }

 	  if(length(geno_vec)>0 ){
 	  	  dis_rd <-c(dis_rd, dis_vec)
		  trio_hap_rd <-c(trio_hap_rd, geno_vec)
 	  	  #print(table(trio_hap_rd))
 	  	  sel <-list()
 	  	  sel <- table(geno_vec[dis_vec<1000000])
 	  	  
 	  	  if(length(sel)>0){
 	  	  		check$prob[cnt] <- 0
 	  	  		check$tol[cnt] <- sum(dis_vec<dis_limit)
 	  	  		check$dis[cnt] <- median(dis_vec[dis_vec<dis_limit])
 	  	  		#print(paste0("Scan ", cnt,":", snv, " ", sum, " ", (p_sum_000 + p_sum_111)/sum))
 	  	  
 	  	  }
 	  }
}



res <-list()
res$dis <- dis_rd
res$hap <- trio_hap_rd
res$check <- check
res$dis2 <- dis_rd2
res$hap2 <- hap_rd2
saveRDS(res, file=out)





