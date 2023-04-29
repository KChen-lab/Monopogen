library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
vega_20 = c(
    '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD',
    '#8C564B', '#E377C2', '#BCBD22', '#17BECF', '#AEC7E8',
    '#FFBB78', '#98DF8A', '#FF9896', '#C5B0D5', '#C49C94',
    '#F7B6D2', '#DBDB8D', '#9EDAE5', '#AD494A', '#8C6D31')

vega_20 <- c(vega_20[4], vega_20[1], vega_20[2], vega_20[3], vega_20[seq(5, length(vega_20),1)])

ref <- read.csv(file=args[1],head=T)
colnames(ref)[1]<-"PopID"
ref <- ref[,seq(1,4,1)]
study <- read.table(file=args[2],head=T,sep="\t")

study_coord <-data.frame("PopID"="study", "indivID"=study[,1], "PC1"=study[,7], "PC2"=study[,8])
dt_plt <- rbind(ref,study_coord)

p <- ggplot(dt_plt, aes(x=PC1, y=PC2, color=PopID)) + geom_point(size=4,alpha=0.1) +  scale_color_manual(values=c(vega_20[1:7],"black")) + theme_bw() + theme( panel.border = element_blank(), panel.grid.major = element_blank(),
panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 16)) + guides(colour=guide_legend(override.aes = list(alpha=0.9)))

p <- p + geom_point(data=dt_plt[dt_plt$PopID=="study",],  col="black",size=4, shape=4)
pdf(file=paste0(args[3],".pdf"),width=4.5, height=3)
print(p)
dev.off()
