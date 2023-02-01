part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsds<-seq(from=0,to=100,by=1)
setwd(part_name)
part<-paste0(part_name,"din/")
setwd(part)
if (!dir.exists("df_RMSD_merged")) {dir.create("df_RMSD_merged")}
df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand))
df_all$center<-NULL
df_all<-unique(df_all)
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("RMSD_merged/",df_all$name[i],".csv"))){
    df_all$receptor[i]<-NA
  }
}
df_all<-df_all%>%filter(!is.na(receptor))
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("df_RMSD_merged/",df_all$name[i],".csv"))){
    df_RMSD_merged<-read.csv(paste0("RMSD_merged/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_RMSD_merged<-df_RMSD_merged%>%mutate(RMSD=round(RMSD,digits = 1))
    df_RMSD_merged<-df_RMSD_merged%>%group_by(RMSD)%>%mutate(number=n())
    df_RMSD_merged<-df_RMSD_merged%>%select(RMSD,number)
    df_RMSD_merged<-unique(df_RMSD_merged)
    write.csv(df_RMSD_merged,paste0("df_RMSD_merged/",df_all$name[i],".csv"),row.names = F)
  }
}

df_RMSD_merged<-read.csv(paste0("df_RMSD_merged/",df_all$name[1],".csv"),stringsAsFactors = F)
for (i in 2:nrow(df_all)) {
  df_RMSD_add<-read.csv(paste0("df_RMSD_merged/",df_all$name[i],".csv"),stringsAsFactors = F)
  df_RMSD_merged<-rbind(df_RMSD_merged,df_RMSD_add)
}
df_RMSD_merged<-df_RMSD_merged%>%group_by(RMSD)%>%mutate(group=sum(number))
df_RMSD_merged<-df_RMSD_merged%>%select(RMSD,group)
df_RMSD_merged<-unique(df_RMSD_merged)
seqv<-seq(from=0,to=max(df_RMSD_merged$RMSD),by=0.1)


seqv<-seqv[!seqv%in%df_RMSD_merged$RMSD]
df_RMSD_add<-data.frame(matrix(ncol=2,nrow = length(seqv)))
colnames(df_RMSD_add)<-colnames(df_RMSD_merged)
df_RMSD_add$RMSD<-seqv
df_RMSD_add$group<-0
df_RMSD_merged<-rbind(df_RMSD_merged,df_RMSD_add)

df_RMSD_merged<-df_RMSD_merged%>%mutate(number=RMSD*10)
df_RMSD_merged<-df_RMSD_merged%>%mutate(numbera=number%/%1)
df_RMSD_merged<-df_RMSD_merged%>%group_by(numbera)%>%mutate(RMSD_new=mean(RMSD))
df_RMSD_merged<-df_RMSD_merged%>%group_by(numbera)%>%mutate(group_new=mean(group))
df_RMSD_merged<-ungroup(df_RMSD_merged)
df_RMSD_merged<-df_RMSD_merged%>%select(RMSD_new,group_new)
df_RMSD_merged<-unique(df_RMSD_merged)
p<-ggplot(data=df_RMSD_merged)+
  labs(x="RMSD, A",y="count")+
  geom_line(aes(x=RMSD_new,y=group_new))+
  geom_point(aes(x=RMSD_new,y=group_new))+
  theme_bw()+
  scale_x_continuous(breaks = v_rmsds,labels = v_rmsds)
ggsave(p,filename = paste0("calibration_RMSD_merge_structure.png"), width = 24, height = 15, units = c("cm"), dpi = 200 )
