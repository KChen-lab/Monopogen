args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(e1071)
library(ggplot2)
options(warn=-1)

mat_gz <- args[1]
outdir <- args[2]
region <- args[3]

### count SNV chuck regions that no germline SNVs included.
SNV_block <-function(summary=NULL){
	block <- 1
	last <- ".|."
	summary$region <- 0
	for(i in seq(1,nrow(summary),1)){
		if(summary[i,10]==".|."){
			if(last==".|." | last=="0|0"){;}
			else{	
				block <- block + 1
			}
			summary$region[i] <- block
		}
		last <- summary[i,10]
	}
	return(summary)
}

### generate positive and negative labels based on position distribution 
SVM_prepare <-function(x=NULL){
	svm<-list()
	region_cnt <- table(x$region)
	filter <- names(region_cnt)[region_cnt>=4]
	neg <- x[(x$region%in% filter & x$genotype==".|."),]
	pos <- x[!x$genotype==".|.",]
	svm$neg <- neg 
	svm$pos <- pos
	svm$test <- x[(!x$region%in%filter) & x$genotype==".|.",]
	return(svm)
}

### SVM trainning on 8 features, some features may not informative (depending on sequencing platforms) 
SVM_train <- function(label=NULL, dir=NULL, region=NULL){

	features <-c("QS", "VDB", "SGB", "RPB", "MQB", "MQSB", "BQB", "MQ0F")
	label$pos <- as.data.frame(label$pos)
	label$pos[label$pos=="None"] <- NA
	# using median values to replace the missing values 
	train_x_pos <- impute(as.matrix(data.matrix(label$pos[,features])), what="median")

	# using the minior value of QS 
	vec <- train_x_pos[,colnames(train_x_pos)=="QS"]
	vec[vec>0.5] <- 1- vec[vec>0.5]
	train_x_pos[,1] <- vec
	
	label$neg <- as.data.frame(label$neg)
	label$neg[label$neg=="None"] <- NA
	train_x_neg <- impute(as.matrix(data.matrix(label$neg[,features])), what="median")
	vec <- train_x_neg[,colnames(train_x_neg)=="QS"]
	vec[vec>0.5] <- 1- vec[vec>0.5]
	train_x_neg[,1] <- vec
	
	label$test <- as.data.frame(label$test)
	label$test[label$test=="None"] <- NA
	test_x <- impute(as.matrix(data.matrix(label$test[,features])), what="median")
	vec <- test_x[,colnames(test_x)=="QS"]
	vec[vec>0.5] <- 1- vec[vec>0.5]
	test_x[,1] <- vec
	train_y_pos <- rep("POS", nrow(train_x_pos))
	train_y_neg <- rep("NEG", nrow(train_x_neg))
	
	# generate feature distribution plot to confirm whether SVM works or not
	pdf(file=paste0(dir,"svm_feature.", region, ".pdf"),width=4,height=5)
	for(i in features){
	  plt_dt <- data.frame("Val"=c(train_x_pos[,i], train_x_neg[,i]),"Label"=c(train_y_pos,train_y_neg))
	  p <- ggplot(aes(fill=Label, y=Val, x=Label),data=plt_dt) + geom_boxplot(position="dodge", alpha=0.5) + 
	    ylab(i) + theme_classic()
	  print(p)
	}
	dev.off()
    
	model <- svm(as.matrix(rbind(train_x_pos, train_x_neg)), 
				 as.factor(c(train_y_pos,train_y_neg)), 
				 probability=TRUE)
	pred <- predict(model, as.matrix(test_x), probability=TRUE)
	prob <- attr(pred, "probabilities")
	prob <-as.data.frame(prob)
	label$test$POS <- prob$POS
	label$test$NEG <- prob$NEG
	return(label)
}


twoloci <- function(mat=NULL, germIndex=NULL, somaticIndex=NULL,dis=NULL){
  
  #germIndex <- germlineIndex
  #somaticIndex <- somaticIndex
  somatic_dis_twoloci <-c()
  somatic_switch_twoloci <-c()
  somatic_posIndex <-c()
  germline_dis_twoloci <-c()
  germline_switch_twoloci <-c()
  
  for(i in seq(1,ncol(mat),1)){
    pos <- intersect(which(!mat[,i]=="0|0"),germIndex)
    if(length(pos)>2){
      somaticIndex_i <- somaticIndex[somaticIndex>min(pos) & somaticIndex<max(pos) & !mat[somaticIndex,i]=="0|0"]
      rd1 <- c()
      rd2 <- c()
      for(j in somaticIndex_i){
            index <- max(which(pos<j))
            lower1 <- pos[index]
            upper1 <- pos[index+1]
            if((dis[j]-dis[lower1]) > (dis[upper1]-dis[j])){
              rd1 <-c(rd1, abs(dis[upper1]-dis[j]))
              rd2 <-c(rd2, as.numeric(mat[j,i]==mat[upper1,i]))
            } else {
              rd1 <-c(rd1, abs(dis[lower1]-dis[j]))
              rd2 <-c(rd2, as.numeric(mat[j,i]==mat[lower1,i]))
            }
      }
      somatic_dis_twoloci <- c(somatic_dis_twoloci, rd1)
      somatic_switch_twoloci <- c(somatic_switch_twoloci, rd2)
      somatic_posIndex <-c(somatic_posIndex, somaticIndex_i)
      ############# extract profiles from germline SNVs only #######################
      lower <- pos[seq(1,length(pos)-1,1)]
      upper <- pos[seq(2,length(pos), 1)]
      germline_dis_twoloci <- c(germline_dis_twoloci, abs(dis[upper]-dis[lower]))
      germline_switch_twoloci <- c(germline_switch_twoloci, as.numeric(mat[lower,i]==mat[upper,i]))
    }
  }
  
  res <- list()
  res$somatic_dis <- somatic_dis_twoloci 
  res$somatic_switch <- somatic_switch_twoloci
  res$somatic_posIndex <- somatic_posIndex
  res$germline_dis <- germline_dis_twoloci 
  res$germline_switch <- germline_switch_twoloci
  return(res)
}

trioloci <- function(mat=NULL, germIndex=NULL, somaticIndex=NULL,dis=NULL){
  
  #germIndex <- germlineIndex
  #somaticIndex <- somaticIndex
  somatic_dis_trioloci <-c()
  somatic_switch_trioloci <-c()
  somatic_posIndex <-c()
  
  germline_dis_trioloci <-c()
  germline_switch_trioloci <-c()

  for(i in seq(1,ncol(mat),1)){
    pos <- intersect(which(!mat[,i]=="0|0"),germIndex)
    if(length(pos)>2){
      somaticIndex_i <- somaticIndex[somaticIndex>min(pos) & somaticIndex<max(pos) & !mat[somaticIndex,i]=="0|0"]
      rd1 <- c()
      rd2 <- c()
      rd3 <- c()
      for(j in somaticIndex_i){
        index <- max(which(pos<j))
        lower1 <- pos[index]
        upper1 <- pos[index+1]
        if(mat[lower1,i]==mat[upper1,i]){
          rd1 <-c(rd1, abs(dis[upper1]-dis[lower1]))
          rd2 <-c(rd2, as.numeric(mat[j,i]==mat[upper1,i]))
          rd3 <-c(rd3, j)
        } 
      }
      somatic_dis_trioloci <- c(somatic_dis_trioloci, rd1)
      somatic_switch_trioloci <- c(somatic_switch_trioloci, rd2)
      somatic_posIndex <-c(somatic_posIndex, rd3)
      
      ############  extract profiles from germline SNVs only #######################
      lower <- pos[seq(1,length(pos)-2,1)]
      middle <- pos[seq(2,length(pos)-1,1)]
      upper <- pos[seq(3,length(pos),1)]
      pos_filter<-which(mat[lower,i]==mat[upper,i])
      if(length(pos_filter)>1){
        upper <- upper[pos_filter]
        middle <- middle[pos_filter]
        lower <- lower[pos_filter]
        germline_dis_trioloci <- c(germline_dis_trioloci, abs(dis[upper]-dis[lower]))
        germline_switch_trioloci <- c(germline_switch_trioloci, as.numeric(mat[lower,i]==mat[middle,i]))
      }
    }
  }
  

  res <- list()
  res$somatic_dis <- somatic_dis_trioloci 
  res$somatic_switch <- somatic_switch_trioloci
  res$somatic_posIndex <- somatic_posIndex
  res$germline_dis <- germline_dis_trioloci 
  res$germline_switch <- germline_switch_trioloci
  return(res)
}

weigthedP <- function(dis=NULL, match=NULL, table=NULL){
  binsize <- c(0, 100,  500, 1000, 2500,  5000, 7500,  10000, 20000, 50000, 100000, 500000, 1000000000000000000)
  rd_p <-c()
  rd_w <-c()
  for(i in seq(2,length(binsize),1)){
    index <-which(dis<binsize[i] & dis>=binsize[i-1])
    if(length(index>0)){
      pval <- sum(match[index])/length(index)
      rd_p <-c(rd_p, pval)
      rd_w <-c(rd_w, (1-table$Prob[i])^2)
    }
  }
  p <- 1-sum(rd_p*rd_w)/(sum(rd_w))
  return(p)
}

calcP <- function(twoloci = NULL, trioloci=NULL, somaticIndex = NULL){
  
  binsize <- c(0, 100,200,300,400, 500, 1000, 2500,  5000, 7500,  10000, 20000, 50000, 100000, 500000, 1000000000000000000)
  table2 <- data.frame("bin"=binsize[seq(2,length(binsize),1)], "tol"=NA, "match"=NA, "Prob"=NA)
  for(i in seq(2,length(binsize),1)){
    pos <-which(twoloci$germline_dis<binsize[i] & twoloci$germline_dis>=binsize[i-1])
    if(length(pos)>100){
      table2$tol[i-1] <- length(pos)
      table2$match[i-1] <- sum(twoloci$germline_switch[pos])
      table2$Prob[i-1] <- 1-table2$match[i-1]/table2$tol[i-1] 
    }
  }
  table3 <- data.frame("bin"=binsize[seq(2,length(binsize),1)], "tol"=NA, "match"=NA, "Prob"=NA)
  for(i in seq(2,length(binsize),1)){
    pos <-which(trioloci$germline_dis<binsize[i] & trioloci$germline_dis>=binsize[i-1])
    if(length(pos)>100){
      table3$tol[i-1] <- length(pos)
      table3$match[i-1] <- sum(trioloci$germline_switch[pos])
      table3$Prob[i-1] <- 1-table3$match[i-1]/table3$tol[i-1] 
    }
  }
  
  # calculate probabilities for putative Somatic SNVs
  # two loci model 
  
  out <- data.frame("SomaticID"=somaticIndex, "p_2"=NA, "p_3"=NA, "p_final"=NA)
  lst <- sort(unique(twoloci$somatic_posIndex))
  for(i in lst){
    pos <- which(twoloci$somatic_posIndex==i)
    dis <- twoloci$somatic_dis[pos]
    match <- twoloci$somatic_switch[pos]
    p2 <- weigthedP(dis=dis, match=match, table=table2)
    out$p_2[out$SomaticID==i] <- p2
  }
  
  lst <- sort(unique(trioloci$somatic_posIndex))
  for(i in lst){
    pos <- which(trioloci$somatic_posIndex==i)
    dis <- trioloci$somatic_dis[pos]
    match <- trioloci$somatic_switch[pos]
    p3 <- weigthedP(dis=dis, match=match, table=table3)
    out$p_3[out$SomaticID==i] <- p3
  }

  for(i in seq(1,nrow(out),1)){
    if(!is.na(out$p_2[i]) &  !is.na(out$p_3[i])){
      out$p_final[i] <- (out$p_2[i] + out$p_3[i])/2
    }
    if(!is.na(out$p_2[i]) &  is.na(out$p_3[i])){
      out$p_final[i] <- out$p_2[i]
    }
    if(is.na(out$p_2[i]) &  !is.na(out$p_3[i])){
      out$p_final[i] <- out$p_3[i]
    }
  }
  res <- list()
  res$table2 <- table2 
  res$table3 <- table3 
  res$somaticLD <- out
  return(res)
}

somaticLD <- function(mat=NULL, svm=NULL,  dir=NULL, region=NULL, min_size=2000){

  mat <- as.data.frame(mat)
  allID <- paste0(mat$V1,mat$V2,mat$V3,mat$V4)
  posID <- paste0(svm$pos$chr,svm$pos$pos,svm$pos$ref,svm$pos$alt)
  testID <- paste0(svm$test$chr,svm$test$pos,svm$test$ref,svm$test$alt)
  mat <-mat[which(allID%in%c(posID, testID)),]
  index <- paste0(mat$V1,mat$V2,mat$V3,mat$V4)
  germlineIndex <- which(index%in%posID)
  somaticIndex <- which(index%in%testID)
  dis <- mat$V2
  # update phasing information 
  #for(i in seq(1,nrow(mat),1)){
  #  if(mat$V10[i]=="1|0"){
  #    tp <- mat[i,]
  #    tp[tp=="1|0"] <- "00|11"
  #    tp[tp=="0|1"] <- "1|0"
  #    tp[tp=="00|11"] <- "0|1"
  #    mat[i,] <- tp
  #  }
  #}
  mat <- mat[,c(seq(19,ncol(mat),1))]
  # initilize the haplotype of somatic SNVs as 0|1
  mat[mat=="0/0"] <- "0|0"
  mat[mat=="0/1"] <- "0|1" 
  mat[mat=="1/0"] <- "1|0"

  n_tol <-nrow(mat)
  if(n_tol > min_size){

    twoloci <- twoloci(mat=mat, germIndex=germlineIndex, somaticIndex = somaticIndex, dis=dis)
    trioloci <- trioloci(mat=mat, germIndex=germlineIndex, somaticIndex = somaticIndex, dis=dis)
    LDrefine <- calcP(twoloci = twoloci, trioloci =  trioloci, somaticIndex = somaticIndex)
    LDrefine$somaticLD$marker <-  index[LDrefine$somaticLD$SomaticID]






    tp <- paste0(svm$test$chr,svm$test$pos,svm$test$ref, svm$test$alt)
    svm$test$p_twoLoci <- NA
    svm$test$p_trioLoci <- NA
    svm$test$p_LDrefine <- NA
    for(i in seq(1,nrow(svm$test),1)){
      pos <- which(LDrefine$somaticLD$marker==tp[i])
      if(length(pos)>0){
        svm$test$p_twoLoci[i] <- LDrefine$somaticLD$p_2[pos]
        svm$test$p_trioLoci[i] <- LDrefine$somaticLD$p_3[pos]
        svm$test$p_LDrefine[i] <- LDrefine$somaticLD$p_final[pos]
      }
    }
    res <- list()
    res$svm <- svm

    # adjust the phasing information based on LD refinement score

    dt <- res$svm$test

    a <- dt$p_twoLoci
    b <- dt$p_trioLoci
    flag_two <- which(a>0.5)
    a[flag_two] <- (1-a[flag_two])
    flag_trio <- which(b>0.5)
    b[flag_trio] <-(1-b[flag_trio])

    # adjust based on their phasing information 
    dt$p_twoLoci <- a
    dt$p_trioLoci <- b
    dt$p_LDrefine <- NA
    for(i in seq(1,nrow(dt),1)){
       if(is.na(a[i])){a[i]<-10}
       if(is.na(b[i])){b[i]<-10}
       dt$p_LDrefine[i] <- (a[i]+b[i])*0.5
       if(dt$p_LDrefine[i]>1){
        dt$p_LDrefine[i] <- NA
       }
    }
    

    baf <- (dt$dep3+dt$dep4)/(dt$dep1+dt$dep2+dt$dep3+dt$dep4)
    res$out <- data.frame("chr"=dt$chr,"pos"=dt$pos,"Ref_allele"=dt$ref,
      "Alt_allele"=dt$alt,"Depth_total"=dt$Dep,"Depth_ref"=dt$dep1+dt$dep2,
      "Depth_alt"=dt$dep3+dt$dep4, "SVM_pos_score"=dt$POS, "LDrefine_twoLoci_score"=dt$p_twoLoci,
      "LDrefine_trioLoci_score"=dt$p_trioLoci, "LDrefine_merged_score"=dt$p_LDrefine,
      "BAF_alt"=baf)
    rownames(res$out) <- paste0(res$out$chr,":",res$out$pos,":", res$out$Ref_allele,":", res$out$Alt_allele)
    res$svm$test <- dt
    res$LDrefine_somatic <- res$out
    res$LDrefine_germline2 <- LDrefine$table2
    res$LDrefine_germline3 <- LDrefine$table3

    table2 <- LDrefine$table2
    table3 <- LDrefine$table3
    table2$model <- "TwoLoci"
    table3$model <- "TrioLoci"
    plt_dt <- rbind(table2, table3)
    plt_dt <- plt_dt[!is.na(plt_dt$tol),]
    plt_dt$bin[plt_dt$bin>5*10^5] <- 50*10^5



    p <- ggplot(plt_dt, aes(x=bin, y=Prob)) +
    geom_point(color="cadetblue",size=6)+
    geom_smooth(color="darksalmon") + scale_x_continuous(trans='log10') +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlab("Physical distance (bp)") +
    ylab("LD refinement score") +
    theme(axis.text=element_text(size=12,angle = 45, vjust = 0.5, hjust=1),
          axis.title=element_text(size=14)) + facet_wrap(~model) + ylim(0,0.65)
    pdf(file=paste0(dir,"LDrefinement_germline.", region, ".pdf"),width=8,height=4)
    print(p)
    dev.off()
    return(res)
  }
  else{
      print(paste0("the No.of germline SNVs in ", region, " are too small (<", min_size, ") to perform LD refinement..."))
      quit(status = 1, save="no")
  }
}





################################### Main function ########################################################
# column information 
# 1: chr  2:pos 3:ref 4:alt 5:dep 6:dep1 7 dep2 8 dep3 9 dep4 
# 10:genotype 11:QS 12:VDB  13:RPB  14:MQB 15:BQB 16:MQSB 17:SGB 18:MQ0F
# 19-end genotype at single cell level 
# dt <- fread("/rsrch1/bcb/kchen_group/jdou/Monopogen/dev/chr20:1000000-2000000.gl.filter.hc.cell.mat.gz")

# outdir <- "/rsrch1/bcb/kchen_group/jdou/Monopogen/dev/"
# region <- "test"

dt <- fread(mat_gz)
dt <- data.frame(dt)
dt$V10[is.na(dt$V10)] <- ".|."
rownames(dt) <- paste0(dt$V1,":",dt$V2,":",dt$V3,":",dt$V4)

cellName <- read.table(file=paste0(outdir,region,".cell.txt"),header=F)

for(j in seq(11,18,1)){
  dt[,j] <- as.numeric(dt[,j])
}
meta <- dt[,1:18]
colnames(meta) <-c("chr","pos","ref","alt","Dep","dep1","dep2","dep3","dep4",
	"genotype","QS","VDB","RPB","MQB","BQB","MQSB","SGB","MQ0F")


mutation_block <- SNV_block(summary=meta)
svm_in <- SVM_prepare(mutation_block)
svm_out <- SVM_train(label =svm_in,dir=outdir, region=region)

final <- somaticLD(mat=dt, svm=svm_out, dir=outdir, region=region)


mut_mat <- dt[rownames(final$out),]
colnames(mut_mat) <-c(colnames(meta),cellName[,1])

write.csv(final$LDrefine_somatic,   paste0(outdir,region,".putativeSNVs.csv"),quote=FALSE,row.names = FALSE)
write.csv(final$LDrefine_germline2, paste0(outdir,region,".germlineTwoLoci_model.csv"),quote=FALSE, row.names = FALSE)
write.csv(final$LDrefine_germline3, paste0(outdir,region,".germlineTrioLoci_model.csv"),quote=FALSE, row.names = FALSE)
saveRDS(mut_mat, file=paste0(outdir,region,".SNV_mat.RDS"))
