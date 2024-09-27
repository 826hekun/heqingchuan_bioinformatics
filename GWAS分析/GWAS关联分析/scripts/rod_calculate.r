#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2021.09.15
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='The reduction of diversity (ROD), defined as: ROD = 1 - π(domesticated)/π(wild)')

parser$add_argument( "-a", "--domesticated", type="character",required=T,
		help="input popA pi result  file for domesticated pop[required]",
		metavar="domesticated")

parser$add_argument( "-b", "--wild", type="character",required=T,
		help="input popB pi result   file for wild pop[required]",
		metavar="wild")

parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
parser$add_argument("-p", "--prefix", type="character", default="ROD",
		help="out file name prefix [default %(default)s]",
		metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}

pop1<-read.table(opt$domesticated,sep="\t",header=T,comment.char="")
pop1<-data.frame(CHROM=pop1$CHROM,BIN_START=pop1$BIN_START, BIN_END=pop1$BIN_END,PIA=pop1$PI)

#pop1$ID=paste(pop1$CHROM,pop1$BIN_START, pop1$BIN_END)
pop2<-read.table(opt$wild,sep="\t",header=T,comment.char="")
pop2<-data.frame(CHROM=pop2$CHROM,BIN_START=pop2$BIN_START, BIN_END=pop2$BIN_END,PIB=pop2$PI)

data<-merge(pop1,pop2,by=c( "CHROM" ,"BIN_START",  "BIN_END"))
#data<-subset(data,data$PIA/data$PIB>1)

data$ROD<-1-(data$PIA/data$PIB)


colnames(data)<-c( "CHROM" ,"BIN_START",  "BIN_END",opt$domesticated,opt$wild,"ROD")
head(data)
write.table(data,
		file=paste0(opt$outdir,"/",opt$prefix,".txt"), 
		sep = "\t",row.names = F,quote =F)
