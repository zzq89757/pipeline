library(dplyr)

args<-commandArgs(TRUE)
#sample=args[1]
#tumor_path=args[2]
querysample=args[1]
T_N_ratio=args[2]
background_path=args[3]
optfinalsample=args[4]
#T_N_ratio=args[3]
#background_path=args[4]

#tumor<-read.table(paste(tumor_path,"/result/",sample,".query.xls",sep=""),sep="\t",header=F)
tumor<-read.table(querysample,sep="\t",header=F)
background<-read.table(paste(background_path,"/background.xls",sep=""),sep="\t",header=T)
background<-background[,c(1:5,21,25,26)]
colnames(background)<-c("index_n","chr_n","pos_n","ref_n","alt_n","AF_background","count_background","total_background")

tumor$index_t<-paste(tumor$V2,tumor$V3,tumor$V4,tumor$V5,sep="_")
colnames(tumor)<-c("sample_t","chr","pos","ref","alt","gene","type","trans","HGVS_c","HGVS_p","AF","exon","DP","SB","qual","ref_F","ref_R","alt_F","alt_R","COSMIC","CLNDN","CLNSIG","avsnp150","ExAC_ALL","ExAC_EAS","1000G_ALL","1000G_EAS","index_t")

tumor_filter<-tumor[!tumor$index_t%in%background$index_n,]
tumor_filter<-tumor_filter[,c(28,1:27)]
tumor_T_N<-merge(tumor,background,by.x="index_t",by.y="index_n")
tumor_T_N$AF<-as.numeric(tumor_T_N$AF)
tumor_T_N$AF_background<-as.numeric(tumor_T_N$AF_background)
tumor_T_N<-tumor_T_N[which(tumor_T_N$AF/tumor_T_N$AF_background>=T_N_ratio),]
tumor_T_N<-tumor_T_N[,c(3:28,33:35)]
tumor_filter<-tumor_filter[,c(3:28)]
tumor_all<-bind_rows(tumor_T_N,tumor_filter)
#write.table(tumor_all,paste(tumor_path,"/result/",sample,"_final.xls",sep=""),sep="\t",row.names=F,col.names=T,quote=F,na="")
write.table(tumor_all,optfinalsample,sep="\t",row.names=F,col.names=T,quote=F,na="")
