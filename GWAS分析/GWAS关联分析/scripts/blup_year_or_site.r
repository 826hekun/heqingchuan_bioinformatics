#!//usr/bin/env Rscript
############################################################################
#北京组学生物科技有限公司
#author huangls
#date 2021.10.28
#version 1.0
#学习R课程：
#R 语言入门与基础绘图：
# https://zzw.xet.tech/s/2G8tHr
#R 语言绘图ggplot2等等：
# https://zzw.xet.tech/s/26edpc

###############################################################################################

library("argparse")
parser <- ArgumentParser(description='trait value to blup for multi year or site data')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input trait result  file[required]",
		metavar="filepath")
#parser$add_argument( "-l", "--lines", type="character",required=T,
#		help="c [required]",
#		metavar="filepath")
#parser$add_argument( "-e", "--env", type="character",required=T,
#		help="envirment lines column name[required]",
#		metavar="filepath")
#parser$add_argument( "-r", "--rep", type="character",required=F,default=NULL,
#		help="rep column name   file[optional default %(default)s]",
#		metavar="filepath")
#

parser$add_argument( "-H", "--height", type="double", default=5,
		help="the height of pic   inches  [default %(default)s]",
		metavar="number")
parser$add_argument("-W", "--width", type="double", default=5,
		help="the width of pic   inches [default %(default)s]",
		metavar="number")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default %(default)s]",
		metavar="path")
#parser$add_argument("-n", "--prefix", type="character", default="trait",
#		help="out file name prefix [default %(default)s]",
#		metavar="prefix")
opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}
###################################################################################





library(lme4)
library(tidyr)
library(ggplot2)
data=read.table(opt$input,head=T,as.is=T,ch=F,sep="\t")

###########################


###########################

if(stringr::str_detect(data[1,1],"<Header")){
	nm<-colnames(data)
	trait.names<-unique(nm[2:length(nm)])
	mytrait.data<-data.frame(matrix(NA,nrow=nrow(data)-1,length(trait.names)))
	colnames(mytrait.data)<-trait.names
	for (trait in trait.names){
		print(trait)
		myselected<-c(nm %in% trait)
		myselected[1]=T
		mydata<-data[,myselected,drop=F]
		colnames(mydata)<-c("<Trait>",rep(trait,ncol(mydata)-1))
		subT<-mydata[2:nrow(mydata),1:ncol(mydata)]
		#print(head(mydata))
		#print(head(subT))
		colnames(subT)<-c("SampleID",mydata[1,2:ncol(mydata)])
		print(head(subT))
		subT=gather(data=subT,key="env",value=!!(trait),2:ncol(subT),na.rm=T)
		print(head(subT))
		#########BLUP ana###########################
		f=as.formula(paste0(trait,"~(1|env)+(1|SampleID)"))
		subT[,trait]=as.double(subT[,trait])
		subT$env=as.factor(subT$env)
		
		blp=lmer(f,data=subT)
		
		blups= ranef(blp)
		
		lines=blups$SampleID+blp@beta
		#print(lines)
		res=data.frame("SampleID"=rownames(lines),blup=lines)
		#print(head(res))

		colnames(res)<-c("<Trait>",trait)
		
		write.table(res,file=paste0(opt$outdir,"/",trait,"_BLUP_value.tsv"),row.names = F,quote = F,sep="\t")
		#hist(lines[,1],col="#0AB3CA",border="white",xlab="BLUP of lines",main="")
		p<-ggplot(res, aes(x=get(trait))) +
				geom_histogram(fill="red")+xlab(trait)+ylab("count")+theme_bw()+ theme(  
						panel.grid=element_blank(), 
						axis.text.x=element_text(colour="black"),
						axis.text.y=element_text(colour="black"),
						panel.border=element_rect(colour = "black"),
						legend.key = element_blank(),
						plot.title = element_text(hjust = 0.5))
		png(filename=paste(opt$outdir,"/",trait,"_BLUP_histogram.png",sep=""), height =opt$height* 300, width = opt$width*300, res = 500, units = "px")
		print(p)
		dev.off()
		
		pdf(file=paste(opt$outdir,"/",trait,"_BLUP_histogram.pdf",sep=""), height = opt$height, width = opt$width)
		print(p)
		dev.off()
		
	}
}else{
	stop("<Header name=env> line not in file ,please check\n")
}









