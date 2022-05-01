part_name <- commandArgs(trailingOnly=TRUE)
#group ligand structures
#part_name<-part_name
library(bio3d)
#library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-1
setwd(part_name)
part<-paste0(part_name,"din/")
setwd(part)
if (dir.exists(paste0("groups/"))) { system(command = paste0("rm -r ",part_name,"din/groups/"))}
if (dir.exists(paste0("groups_fin/"))) {system(command = paste0("rm -r ",part_name,"din/groups_fin/"))}
if (dir.exists(paste0("str_fin/"))) {system(command = paste0("rm -r ",part_name,"din/str_fin/"))}

if (!dir.exists("groups")) {dir.create("groups")}

if (!dir.exists("RMSD_analysis/")){dir.create("RMSD_analysis/")}
df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
i<-1
v_RMSD_analysis<-list.files("RMSD_analysis/")
i<-1
if(length(v_RMSD_analysis)>0){
  a<-c()
  for (i in 1:length(v_RMSD_analysis)) {
    b<-strsplit(v_RMSD_analysis[i],split = ".",fixed = T)[[1]][1]
    a<-c(a,b)
  }
  v_RMSD_analysis<-a
  df_all<-df_all[!df_all$name%in%v_RMSD_analysis,]
}
if(nrow(df_all)>0){
  for (i in 1:nrow(df_all)) {
    models<-list.files(paste0("pdb_second/",df_all$name[i]))
    if(length(models)>1){
      df_RMSD<-data.frame(matrix(ncol = 2,nrow=length(models)))
      colnames(df_RMSD)<-c("models","RMSD")
      df_RMSD$models<-models
      df_RMSD_all<-full_join(df_RMSD,df_RMSD,by="RMSD")
      df_RMSD_all<-df_RMSD_all%>%filter(models.x!=models.y)
      for (j in 1:nrow(df_RMSD_all)) {
        pdb_1<-read.pdb(paste0("pdb_second/",df_all$name[i],"/",df_RMSD_all$models.x[j]))
        pdb_2<-read.pdb(paste0("pdb_second/",df_all$name[i],"/",df_RMSD_all$models.y[j]))
        df_RMSD_all$RMSD[j]<-rmsd(pdb_1,pdb_2)
      }
      write.csv(df_RMSD_all,paste0("RMSD_analysis/",df_all$name[i],".csv"),row.names = F)
    }
  }
}
