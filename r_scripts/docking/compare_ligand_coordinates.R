part_name <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)

df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%select(receptor,ligand)
df_all<-unique(df_all)
part_name<-paste0(part_name,"din/")
setwd(part_name)

df_merge_surf<-read.csv(paste0(part_name,"df_merge_structure_log.csv"),stringsAsFactors = F)
df_merge_center<-read.csv(paste0(part_name,"df_merge_structure_log_center.csv"),stringsAsFactors = F)

df_merge_surf<-df_merge_surf%>%select(name.x,receptor,ligand)
df_merge_surf<-unique(df_merge_surf)

df_merge_center<-df_merge_center%>%select(name.x,receptor,ligand)
df_merge_center<-unique(df_merge_center)
df_merge<-rbind(df_merge_center,df_merge_surf)

#df_merge<-read.csv(file = paste0(part_name,"RMSD_comparition.csv"),stringsAsFactors = F)
#df_merge<-df_merge%>%filter(RMSD<5)
df_merge<-df_merge%>%mutate(x=NA)
df_merge<-df_merge%>%mutate(y=NA)
df_merge<-df_merge%>%mutate(z=NA)

for (i in 1:nrow(df_merge)) {
  pdb<-read.pdb(paste0(part_name,"str_fin/",df_merge$name.x[i]))
  df_pdb<-pdb$atom
  df_merge$x[i]<-median(df_pdb$x)
  df_merge$y[i]<-median(df_pdb$y)
  df_merge$z[i]<-median(df_pdb$z)
}
df_merge<-left_join(df_merge,df_merge,by=c("receptor","ligand"),relationship = "many-to-many")
df_merge<-df_merge%>%filter(abs(z.x-z.y)<8)
df_merge<-df_merge%>%mutate(RMSD=NA)
for (i in 1:nrow(df_merge)) {
  pdb_1<-read.pdb(paste0("str_fin/",df_merge$name.x.y[i]))
  pdb_2<-read.pdb(paste0("str_fin/",df_merge$name.x.x[i]))
  df_merge$RMSD[i]<-rmsd(pdb_1,pdb_2)
}
#df_merge<-df_merge%>%filter(RMSD<2.5)
write.csv(df_merge,file = paste0(part_name,"RMSD_comparition.csv"),row.names = F)
