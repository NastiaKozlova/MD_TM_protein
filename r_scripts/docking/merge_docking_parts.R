part_start <- commandArgs(trailingOnly=TRUE)
#group ligand structures
#part_start<-part_name
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4
setwd(part_start)
part<-paste0(part_start,"din/")
setwd(part)

if (!dir.exists("groups_merged")) {dir.create("groups_merged")}

df_all<-read.csv(paste0(part_start,"df_all.csv"),stringsAsFactors = F)

if (!dir.exists("RMSD_analysis/")){dir.create("RMSD_analysis/")}
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))

v_protein_name<-unique(df_all$receptor)
v_structure_RMSD<-list.files(paste0("str_fin/"))
df_structure_RMSD<-data.frame(matrix(ncol=4,nrow = length(v_structure_RMSD)))
colnames(df_structure_RMSD)<-c("name","receptor","ligand","RMSD")
df_structure_RMSD$name<-v_structure_RMSD
j<-1
for (j in 1:nrow(df_structure_RMSD)) {
  df_structure_RMSD$receptor[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][1]
  df_structure_RMSD$ligand[j]<-strsplit(x = df_structure_RMSD$name[j],split = "_")[[1]][2]
}
df_structure_RMSD_analysis<-left_join(df_structure_RMSD,df_structure_RMSD,by=c("receptor","ligand","RMSD"))
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(name.x!=name.y)
j<-1
for (j in 1:nrow(df_structure_RMSD_analysis)) {
  pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.x[j]))
  pdb_2<-read.pdb(paste0("str_fin/",df_structure_RMSD_analysis$name.y[j]))
  
  df_structure_RMSD_analysis$RMSD[j]<-rmsd(pdb_1,pdb_2)
}

df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(RMSD<4)
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(models.x)%>%mutate(number=n())
df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)                                             
df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(number>5)
if (nrow(df_structure_RMSD_analysis)>0){
  df_RMSD<-df_structure_RMSD_analysis%>%select(models.x,number)
  df_RMSD<-unique(df_RMSD)
  df_RMSD<-df_RMSD%>%arrange(desc(number))
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis
  for (j in 1:nrow(df_RMSD)) {
    if (!is.na(df_RMSD$models.x[j])) {
      df_structure_RMSD_analysis_test<-df_structure_RMSD_analysis%>%filter(models.x==df_RMSD$models.x[j])
      if (nrow(df_structure_RMSD_analysis_test)>5) {
        df_structure_RMSD_analysis_test<-df_structure_RMSD_analysis_test%>%mutate(grop_number=j)
        write.csv(df_structure_RMSD_analysis_test,paste0("groups_merged/",df_all$name[i],"/grop_",j,".csv"),row.names = F) 
      }
    }
    
    df_structure_RMSD_analysis$models.x[df_structure_RMSD_analysis$models.x%in%c(df_structure_RMSD_analysis_test$models.y,df_structure_RMSD_analysis_test$models.x)]<-NA
    df_structure_RMSD_analysis$models.x[df_structure_RMSD_analysis$models.y%in%c(df_structure_RMSD_analysis_test$models.y,df_structure_RMSD_analysis_test$models.x)]<-NA
    df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(!is.na(models.x))
    
    df_RMSD$models.x[df_RMSD$models.x%in%df_structure_RMSD_analysis_test$models.y]<-NA
    df_RMSD$models.x[df_RMSD$models.x%in%df_structure_RMSD_analysis_test$models.x]<-NA
  }
}