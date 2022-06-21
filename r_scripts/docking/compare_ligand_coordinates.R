part_analysis <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)

df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%select(receptor,ligand)
df_all<-unique(df_all)
part_name<-paste0(part_analysis,"din/")
setwd(part_name)

df_merge_surf<-read.csv(paste0(part_name,"df_merge_structure_log.csv"),stringsAsFactors = F)
df_merge_center<-read.csv(paste0(part_name,"df_merge_structure_log_center.csv"),stringsAsFactors = F)

df_merge_surf<-df_merge_surf%>%select(name.x,receptor,ligand)
df_merge_surf<-unique(df_merge_surf)

df_merge_center<-df_merge_center%>%select(name.x,receptor,ligand)
df_merge_center<-unique(df_merge_center)

df_merge<-left_join(df_merge_surf,df_merge_center,by=c("receptor","ligand"))
df_merge<-df_merge%>%mutate(RMSD=NA)
for (i in 1:nrow(df_merge)) {
  pdb_1<-read.pdb(paste0("structure_merged_center/",df_merge$name.x.y[i]))
  pdb_2<-read.pdb(paste0("structure_merged/",df_merge$name.x.x[i]))
  df_merge$RMSD[i]<-rmsd(pdb_1,pdb_2)
}
#df_merge<-df_merge%>%filter(RMSD<2.5)
write.csv(df_merge,file = paste0(part_name,"RMSD_comparition.csv"),row.names = F)