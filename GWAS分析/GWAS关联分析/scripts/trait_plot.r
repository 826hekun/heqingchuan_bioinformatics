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
parser <- ArgumentParser(description='trait or phenotype  plot')

parser$add_argument( "-i", "--input", type="character",required=T,
		help="input trait result  file[required]",
		metavar="filepath")
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

library(ggplot2)
library(reshape2)
library("gplots")
library(pheatmap)
library(lme4)
###step 1  plot
od=opt$outdir
data=read.table(opt$input,head=T,as.is=T,ch=F,sep="\t")
#head(data)
color_map<-colorRampPalette(c("blue", "white", "red"))(100)
if(stringr::str_detect(data[1,1],"<Header")){
	  #path=file.path(od,"split")
	  #if(!file.exists(path)){dir.create(path,recursive =T)}
	  nm<-colnames(data)
	  trait.names<-unique(nm[2:length(nm)])
	  mytrait.data<-data.frame(matrix(NA,nrow=nrow(data)-1,length(trait.names)))
	  colnames(mytrait.data)<-trait.names
	  for (trait in trait.names){
		  myselected<-c(nm %in% trait)
		  myselected[1]=T
		  mydata<-data[,myselected,]
		  ##print(head(mydata))
		  
		  colnames(mydata)<-c("<Trait>",rep(trait,ncol(mydata)-1))
		  subT<-mydata[2:nrow(mydata),2:ncol(mydata)]
		  
		  #print(head(t(subT)))
		  d<-apply(as.data.frame(lapply(subT,as.numeric)),1,mean,na.rm=T)
		  mytrait.data[trait]<-d
		  #d=d[!is.na(d)]
		  p<-ggplot(data=as.data.frame(d), aes(d)) +geom_histogram(fill="red")+
				  xlab(trait)+ylab("count")+theme_bw()+ theme(  
						  panel.grid=element_blank(), 
						  axis.text.x=element_text(colour="black"),
						  axis.text.y=element_text(colour="black"),
						  panel.border=element_rect(colour = "black"),
						  legend.key = element_blank(),
						  plot.title = element_text(hjust = 0.5))
		  png(filename=paste(od,"/",trait,"_histogram.png",sep=""), height = opt$height*300, width = opt$width*300, res = 500, units = "px")
		  print(p)
		  dev.off()
		  pdf(file=paste(od,"/",trait,"_histogram.pdf",sep=""), height = opt$height, width = opt$width)
		  print(p)
		  dev.off()
		  #print(quantile(d))
		  write.table(mydata,file=paste(od,"/",trait,".txt",sep=""),row=F,sep="\t",quote=F)
		  OutVals <- boxplot(d, plot=FALSE)$out
		  rmOuters<-c(F, d %in% OutVals)
		  mydata[rmOuters,2:ncol(mydata)]<-NA
		  write.table(mydata,file=paste(od,"/",trait,".rm_outers.txt",sep=""),row=F,sep="\t",quote=F)
	  }
		  #print(head(data))
		  data_scale=mytrait.data
		  #print("aa")
		  data_scale=apply(data_scale,2,scale)
		  #print("bb")
		  d=melt(data_scale)
		  d=data.frame(d)
		  colnames(d)[2:3]=c("x","y")
		  d$x=as.factor(d$x)
		  
#		  p=ggplot(data=d, aes(x=x, y=y, fill=x))+
#				 # geom_violin(trim=FALSE) +
#				  geom_boxplot(notch = FALSE,outlier.size = -1, color="black",lwd=0.5, alpha = 0.7,show.legend = F)+
#				  theme_bw()+
#				  theme( 
#						  #axis.text.x =element_text(angle=45,vjust=0.5,hjust=0.5),
#						  panel.grid=element_blank(), 
#						  panel.grid.minor = element_blank(), 
#						  panel.grid.major = element_blank(), 
#						  legend.position="none" ,
#						  axis.text.x=element_text(angle=45,hjust=1)
#				  )+
#				  labs( 
#						  x='Trait', 
#						  y="Normalization  Value" 
#				  ) 
#		  
#		  png(filename=paste(od,"/","trait_boxplot.png",sep=""), height = opt$height*300, width = opt$width*300, res = 500, units = "px")
#		  print(p)
#		  dev.off()
#		  
#		  pdf(file=paste(od,"/","trait_boxplot.pdf",sep=""), height = opt$height, width = opt$width)
#		  print(p)
#		  dev.off()
		  #print(data_scale)
		  if(ncol(data_scale)>1){
			  ##step2  cor cluster
			  cor=cor(data_scale,use="pairwise.complete.obs")
			  png(filename=paste(od,"/trait_cor_cluster.png",sep=""),width=8*300,height=8*300,units="px",res=300)
			  pheatmap(cor,color=color_map,
					scale="none",
					cluster_row=T,
					cluster_col=T,
					border=FALSE
					  )
			  dev.off()
			  pdf(file=paste(od,"/trait_cor_cluster.pdf",sep=""),width=8,height=8)
			  pheatmap(cor,color=color_map,
					  scale="none",
					  cluster_row=T,
					  cluster_col=T,
					  border=FALSE
			  )
			  dev.off()
			  cor=data.frame(Trait=rownames(cor),cor)
			  
			  write.table(cor,file=paste(od,"/trait_cor.txt",sep=""),sep="\t",row=F,quote=F)
		  }
  }else{
	  
	  
	  print(head(data))
	  for(i in 2:ncol(data)){
		  x=data[,i]
		  x[x==-999]=NA
		  data[,i]=x
	  }
	  path=file.path(od,"split")
	  if(!file.exists(path)){dir.create(path,recursive =T)}
	  nm=colnames(data)
	  for(i in 2:ncol(data)){
		  d=data[,c(1,i)]
		  d=d[!is.na(d[,2]),]
		#print(d[,2])
		  p<-ggplot(d, aes(d[,2])) +
				  geom_histogram(fill="red")+xlab(nm[i])+ylab("count")+theme_bw()+ theme(  
						  panel.grid=element_blank(), 
						  axis.text.x=element_text(colour="black"),
						  axis.text.y=element_text(colour="black"),
						  panel.border=element_rect(colour = "black"),
						  legend.key = element_blank(),
						  plot.title = element_text(hjust = 0.5))
		  png(filename=paste(path,"/",nm[i],"_histogram.png",sep=""), height =opt$height* 300, width = opt$width*300, res = 500, units = "px")
		  print(p)
		  dev.off()
		  
		  pdf(file=paste(path,"/",nm[i],"_histogram.pdf",sep=""), height = opt$height, width = opt$width)
		  print(p)
		  dev.off()
		  print(quantile(d[,2]))
		  write.table(d,file=paste(path,"/",nm[i],".txt",sep=""),row=F,sep="\t",quote=F)
		  OutVals <- boxplot(d[,2], plot=FALSE)$out
		  write.table(d[! d[,2] %in% OutVals,],file=paste(path,"/",nm[i],".rm_outers.txt",sep=""),row=F,sep="\t",quote=F)
	  }
	  #print(head(data))
	  data_scale=data[,-1,drop=F]
	  #print("dd")
	  data_scale=apply(data_scale,2,scale)
	  #print("cc")
	  d=melt(data_scale)
	  d=data.frame(d)
	  colnames(d)[2:3]=c("x","y")
	  d$x=as.factor(d$x)
#	  p=ggplot(data=d, aes(x=x, y=y, fill=x))+
#			  #geom_violin(trim=FALSE) +
#			  geom_boxplot(notch = FALSE,outlier.size = 0.1, color="black",lwd=0.5, alpha = 0.7,show.legend = F)+
#			  theme_bw()+
#			  theme( 
#					  #axis.text.x =element_text(angle=45,vjust=0.5,hjust=0.5),
#					  panel.grid=element_blank(), 
#					  panel.grid.minor = element_blank(), 
#					  panel.grid.major = element_blank(), 
#					  legend.position="none" ,
#					  axis.text.x=element_text(angle=45,hjust=1)
#			  )+labs( 
#					  x='Trait', 
#					  y="Normalization  Value" 
#			  ) 
#	  png(filename=paste(od,"/","trait_boxplot.png",sep=""), height = opt$height*300, width = opt$width*300, res = 500, units = "px")
#	  print(p)
#	  dev.off()
#	  
#	  pdf(file=paste(od,"/","trait_boxplot.pdf",sep=""), height = opt$height, width = opt$width)
#	  print(p)
#	  dev.off()
	  if(ncol(data_scale)>1){
		  
		  
		  ##step2  cor cluster
		  cor=cor(data_scale,use="pairwise.complete.obs")
		  png(filename=paste(od,"/trait_cor_cluster.png",sep=""),width=8*300,height=8*300,units="px",res=300)
			pheatmap(cor,color=color_map,
					scale="none",
					cluster_row=T,
					cluster_col=T,
					border=FALSE
			)
		  dev.off()
		  pdf(file=paste(od,"/trait_cor_cluster.pdf",sep=""),width=8,height=8)
		  pheatmap(cor,color=color_map,
				  scale="none",
				  cluster_row=T,
				  cluster_col=T,
				  border=FALSE
		  )
		  dev.off()
		  cor=data.frame(Trait=rownames(cor),cor)
		  write.table(cor,file=paste(od,"/trait_cor.txt",sep=""),sep="\t",row=F,quote=F)
		  
		  
	  }
	  
  }
 




