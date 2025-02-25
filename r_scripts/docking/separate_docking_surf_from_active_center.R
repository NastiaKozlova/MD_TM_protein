part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
#v_rmsd<-3.5

print(Sys.time())

setwd(part_analysis)
df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
#df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_all<-df_all%>%mutate(x=NA)
df_all<-df_all%>%mutate(y=NA)
df_all<-df_all%>%mutate(z=NA)
i<-1
for (i in 1:nrow(df_all)) {
  a<-strsplit(df_all$center[i],split = "_")[[1]]
  df_all$x[i]<-as.numeric(a[3])
  df_all$y[i]<-as.numeric(a[5])
  df_all$z[i]<-as.numeric(a[7])
}
df_all<-df_all%>%mutate(surf="surf")
df_all$surf[is.na(df_all$x)]<-"center"
df_all<-df_all%>%select(name,receptor,ligand,center,surf,x,y,z)
write.csv(df_all,"df_separate.csv",row.names=F)