
#background variants statistics 
#input:file contain all background samples with titles below:
#sample	CHROM   POS     REF     ALT     ANN[*].GENE     ANN[*].EFFECT   ANN[*].FEATUREID        ANN[*].HGVS_C   ANN[*].HGVS_P   AF      ANN[*].RANK     DP      SB      QUAL    DP4[0]  DP4[1]  DP4[2]  DP4[3]  cosmic92_coding   CLNDN   CLNSIG  avsnp150        EXAC_all  EXAC_eas  ALL.sites.2015_08       EAS.sites.2015_08


args<-commandArgs(TRUE)
input<-args[1]
output<-args[2]

library(dplyr)
all<-read.table(input,sep="\t",header=T,check.names=F)
all$index<-paste(all$CHROM,all$POS,all$REF,all$ALT,sep="_")
all_freq<-all%>%group_by(index)%>%summarize(count=n())
all_freq_count<-all_freq%>%group_by(count)%>%summarize(number=n())

all_AF_min<-all%>%group_by(index)%>%summarise(min(AF))
all_AF_max<-all%>%group_by(index)%>%summarise(max(AF))
all_AF_median<-all%>%group_by(index)%>%summarise(median(AF))
all_DP_median<-all%>%group_by(index)%>%summarise(median(DP))
all_QUAL_median<-all%>%group_by(index)%>%summarise(median(QUAL))
all_unique<-all[!duplicated(all$index),]

final<-merge(all_unique[,c(28,2:10,12,20:27)],all_AF_min,by.x="index",by.y="index")
final<-merge(final,all_AF_median,by.x="index",by.y="index")
final<-merge(final,all_AF_max,by.x="index",by.y="index")
final<-merge(final,all_DP_median,by.x="index",by.y="index")
final<-merge(final,all_QUAL_median,by.x="index",by.y="index")
final<-merge(final,all_freq,by.x="index",by.y="index")
final$total<-length(all_freq_count$count)
write.table(final,output,sep="\t",quote=F,row.names=F)
