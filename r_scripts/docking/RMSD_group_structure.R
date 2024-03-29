part_name <- commandArgs(trailingOnly=TRUE)
#group ligand structures
v_test<-strsplit(part_name,split = ",")[[1]]
part_start<-v_test[1]
lig<-v_test[2]
lig<-lig[!is.na(lig)]
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
#v_rmsd<-4
setwd(part_start)
part<-paste0(part_start,"din/")
setwd(part)

if (!dir.exists("groups")) {dir.create("groups")}

df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)
if (!dir.exists("RMSD_analysis/")){dir.create("RMSD_analysis/")}
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_all<-df_all%>%mutate(number=0)
i<-1
print(Sys.time())
for (i in 1:nrow(df_all)) {
  models<-list.files(paste0("pdb_second/",df_all$name[i]))
  df_all$number[i]<-length(models)
}
df_all<-df_all%>%filter(number>0)

print(Sys.time())
v_str<-list.files(paste0("RMSD_analysis/"))
a<-c()
i<-1
if(length(v_str)>0){
  for (i in 1:length(v_str)) {
    b<-strsplit(v_str[i],split = ".",fixed = T)[[1]][1]
    a<-c(a,b)
  }
  v_str<-a
  df_all<-df_all[!df_all$name%in%v_str,]
}
if(length(lig)==1){
  df_all<-df_all%>%filter(ligand==lig)
}

if(nrow(df_all)>0){
  for (i in 1:nrow(df_all)) {
    if(!file.exists(paste0("RMSD_analysis/",df_all$name[i],".csv"))){
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
}
