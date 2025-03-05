part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
part_TEMP<-strsplit(part_analysis,split = ",")[[1]]
#part_start<-part_TEMP[1]
v_rmsd_start<-as.numeric(part_TEMP[2])
part_analysis<-part_TEMP[1]
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)

setwd(part_analysis)
print("kek")
df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))

part<-paste0(part_analysis,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged_center/"))) {system(command = paste0("rm -r ",part,"fin_merged_center/"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged_center/"))) {system(command = paste0("rm -r ",part,"structure_merged_center/"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged_center/"))) {system(command = paste0("rm -r ",part,"groups_merged_center/"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged_center/")) {dir.create("RMSD_merged_center/")}
if (!dir.exists("groups_merged_center/")) {dir.create("groups_merged_center/")}
if (!dir.exists("structure_merged_center/")) {dir.create("structure_merged_center/")}
if (!dir.exists("fin_merged_center/")) {dir.create("fin_merged_center/")}

df_analysis<-df_all%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
q<-1
for (q in 1:nrow(df_analysis)) {
  if(!file.exists(paste0("RMSD_merged_center/",df_analysis$receptor_ligand[q],".csv"))){
    df_analysis$receptor_ligand[q]<-NA
  }
}
df_analysis<-df_analysis%>%filter(!is.na(receptor_ligand))
q<-1

for (q in 1:nrow(df_analysis)) {
  
  df_structure_RMSD_analysis<-read.csv(paste0("RMSD_merged_center/",df_analysis$receptor_ligand[q],".csv"),stringsAsFactors = F)
  
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(name.x=paste0(ligand_center.x,"_",models.x.x))
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(name.y=paste0(ligand_center.y,"_",models.x.y))
  #fix RMSD threshold
  v_rmsd_temp<-quantile(df_structure_RMSD_analysis$RMSD,probs=0.025)
  if(v_rmsd_start<v_rmsd_temp){v_rmsd<-v_rmsd_temp}else{v_rmsd<-v_rmsd_start}
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%filter(RMSD<v_rmsd)
  df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
  df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(number=n())
  df_structure_RMSD_analysis<-ungroup(df_structure_RMSD_analysis)       
  if (!dir.exists(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))) {
    dir.create(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))
  }
  
  df_structure_RMSD_complex<-df_structure_RMSD_analysis
  if (nrow(df_structure_RMSD_complex)>0){
    
    df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
    df_RMSD<-unique(df_RMSD)
    df_RMSD<-df_RMSD%>%arrange(desc(number))
    v_structure<-c()#unique(df_structure_RMSD_analysis$name.x,df_structure_RMSD_analysis$name.y)
    df_structure_RMSD_complex<-df_structure_RMSD_complex
    temp<-nrow(df_RMSD)
    for (j in 1:nrow(df_RMSD)) {
      if (!is.na(df_RMSD$name.x[j])) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex%>%filter(name.x==df_RMSD$name.x[j])
        if (nrow(df_structure_RMSD_complex_test)>1) {
          df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j)
          v_structure_TEMP<-unique(c(df_structure_RMSD_complex_test$name.x,df_structure_RMSD_complex_test$name.y))
          v_structure<-c(v_structure,v_structure_TEMP)
          write.csv(df_structure_RMSD_complex_test,paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/grop_",j,".csv"),row.names = F) 
        }
      }
      
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.x%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex$name.x[df_structure_RMSD_complex$name.y%in%c(df_structure_RMSD_complex_test$name.y,df_structure_RMSD_complex_test$name.x)]<-NA
      df_structure_RMSD_complex<-df_structure_RMSD_complex%>%filter(!is.na(name.x))
      
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.y]<-NA
      df_RMSD$name.x[df_RMSD$name.x%in%df_structure_RMSD_complex_test$name.x]<-NA
    }
    
    
    df_structure_RMSD_complex<-df_structure_RMSD_analysis
    df_structure_RMSD_complex<-df_structure_RMSD_complex[!df_structure_RMSD_complex$name.x%in%v_structure,]
    df_structure_RMSD_complex<-df_structure_RMSD_complex[!df_structure_RMSD_complex$name.y%in%v_structure,]
    if(nrow(df_structure_RMSD_complex)>0){
      df_RMSD<-df_structure_RMSD_complex%>%select(name.x,number)
      df_RMSD<-unique(df_RMSD)
      df_RMSD<-df_RMSD%>%arrange(desc(number))
      for (j in (1:nrow(df_RMSD))) {
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex[j,]
        df_structure_RMSD_complex_test<-df_structure_RMSD_complex_test%>%mutate(grop_number=j+temp)
        write.csv(df_structure_RMSD_complex_test,paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/grop_",j+temp,".csv"),row.names = F) 
      }
    }
  }
  v_structure<-list.files(paste0("groups_merged_center/",df_analysis$receptor_ligand[q]))
  df_structure_RMSD_analysis_start<-read.csv(paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/",v_structure[1]))
  df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(RMSD<0)
  for (j in 1:length(v_structure)) {
    df_structure_RMSD_analysis<-read.csv(paste0("groups_merged_center/",df_analysis$receptor_ligand[q],"/",v_structure[j]))
    df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
  }
  
  write.csv(df_structure_RMSD_analysis_start,paste0("fin_merged_center/",df_analysis$receptor_ligand[q],".csv"),row.names = F)
}
df_structure_RMSD_analysis_start<-read.csv(paste0("fin_merged_center/",df_analysis$receptor_ligand[1],".csv"),stringsAsFactors = F)
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%group_by(name.x)%>%mutate(size_of_group=n())
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(is.na(name.x))
j<-1
if(nrow(df_analysis)>0){
  for (j in 1:nrow(df_analysis)) {
    if(file.exists(paste0("fin_merged_center/",df_analysis$receptor_ligand[j],".csv"))){
      df_structure_RMSD_analysis<-read.csv(paste0("fin_merged_center/",df_analysis$receptor_ligand[j],".csv"),stringsAsFactors = F)
      df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
      df_structure_RMSD_analysis<-unique(df_structure_RMSD_analysis)
      df_structure_RMSD_analysis_add<-df_structure_RMSD_analysis[1,]
      df_structure_RMSD_analysis_add<-df_structure_RMSD_analysis_add%>%mutate(RMSD=0)
      df_structure_RMSD_analysis_add<-df_structure_RMSD_analysis_add%>%mutate(name.y=name.x)
      df_structure_RMSD_analysis<-rbind(df_structure_RMSD_analysis,df_structure_RMSD_analysis_add)
      df_structure_RMSD_analysis<-df_structure_RMSD_analysis%>%group_by(name.x)%>%mutate(size_of_group=n())
      df_structure_RMSD_analysis_start<-rbind(df_structure_RMSD_analysis_start,df_structure_RMSD_analysis)
    }
  }
}
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%filter(!is.na(name.x))
df_structure_RMSD_analysis_start<-ungroup(df_structure_RMSD_analysis_start)
df_structure_RMSD_analysis_start<-df_structure_RMSD_analysis_start%>%select(name.x, models.x.x, ligand_center.x, center.x, receptor, 
                                                                            name.y, models.x.y, ligand_center.y, center.y, RMSD, 
                                                                            ligand, surf, number, grop_number, size_of_group)
colnames(df_structure_RMSD_analysis_start)
df_log<-read.csv("log_fin.csv",stringsAsFactors = F)
df_log<-df_log%>%select(models.x,models.y,grop_number,ligand,affinity,center,receptor,ligand)

df_log<-df_log%>%mutate(name=paste0(receptor,"_", ligand,"_",center,"_",models.x))
df_structure_RMSD<-left_join(df_structure_RMSD_analysis_start,df_log,by=c("name.y" ="name","ligand","receptor"),
                             relationship = "many-to-many")
write.csv(df_structure_RMSD,"df_merge_structure_log.csv",row.names = F)
setwd(part)
df_structure_RMSD<-read.csv("df_merge_structure_log.csv",stringsAsFactors = F)

#df_structure_RMSD<-df_structure_RMSD%>%filter(RMSD==0)
df_structure_RMSD<-df_structure_RMSD%>%select(name.x,receptor,ligand)
df_structure_RMSD<-unique(df_structure_RMSD)

for (j in 1:nrow(df_structure_RMSD)) {
  pdb_1<-read.pdb(paste0("str_fin/",df_structure_RMSD$name.x[j]))  
  write.pdb(pdb_1,paste0("structure_merged_center/",df_structure_RMSD$name.x[j]))
}
