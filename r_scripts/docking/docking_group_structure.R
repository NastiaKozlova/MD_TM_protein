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

df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
#sort to grops
i<-1
for (i in 1:nrow(df_all)) {
  
  if(file.exists(paste0("RMSD_analysis/",df_all$name[i],".csv"))){
    if (!dir.exists(paste0("groups/",df_all$name[i]))) {dir.create(paste0("groups/",df_all$name[i]))}
    df_RMSD_all<-read.csv(paste0("RMSD_analysis/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_RMSD_all<-df_RMSD_all%>%filter(RMSD<v_rmsd)
    df_RMSD_all<-df_RMSD_all%>%group_by(models.x)%>%mutate(number=n())
    df_RMSD_all<-ungroup(df_RMSD_all)                                             
    df_RMSD_all<-df_RMSD_all%>%filter(number>5)
    if (nrow(df_RMSD_all)>0){
      df_RMSD<-df_RMSD_all%>%select(models.x,number)
      df_RMSD<-unique(df_RMSD)
      df_RMSD<-df_RMSD%>%arrange(desc(number))
      df_RMSD_all<-df_RMSD_all
      for (j in 1:nrow(df_RMSD)) {
        if (!is.na(df_RMSD$models.x[j])) {
          df_RMSD_all_test<-df_RMSD_all%>%filter(models.x==df_RMSD$models.x[j])
          if (nrow(df_RMSD_all_test)>5) {
            df_RMSD_all_test<-df_RMSD_all_test%>%mutate(grop_number=j)
            write.csv(df_RMSD_all_test,paste0("groups/",df_all$name[i],"/grop_",j,".csv"),row.names = F) 
          }
        }
        
        df_RMSD_all$models.x[df_RMSD_all$models.x%in%c(df_RMSD_all_test$models.y,df_RMSD_all_test$models.x)]<-NA
        df_RMSD_all$models.x[df_RMSD_all$models.y%in%c(df_RMSD_all_test$models.y,df_RMSD_all_test$models.x)]<-NA
        df_RMSD_all<-df_RMSD_all%>%filter(!is.na(models.x))
        
        df_RMSD$models.x[df_RMSD$models.x%in%df_RMSD_all_test$models.y]<-NA
        df_RMSD$models.x[df_RMSD$models.x%in%df_RMSD_all_test$models.x]<-NA
      }
    }
  }
}


#combine all groups logs files
i<-1
j<-1
if (!dir.exists("groups_fin")) {dir.create("groups_fin")}
for (i in 1:nrow(df_all)) {
  v_groups<-list.files(paste0("groups/",df_all$name[i]))
  if(length(v_groups)>0){
    df_groups<-read.csv(paste0("groups/",df_all$name[i],"/",v_groups[1]),stringsAsFactors = F)
    df_groups_TEMP<-df_groups%>%mutate(models.y=models.x)
    df_groups_TEMP<-df_groups_TEMP%>%mutate(RMSD=0)
    df_groups_TEMP<-unique(df_groups_TEMP)
    df_groups<-rbind(df_groups,df_groups_TEMP)
    
    df_groups<-df_groups%>%mutate(group=v_groups[1])
    df_groups<-df_groups%>%mutate(ligand_center=df_all$name[i])
    
    if (length(v_groups)>1) {
      for (j in 2:length(v_groups)) {
        df_groups_add<-read.csv(paste0("groups/",df_all$name[i],"/",v_groups[j]))
        df_groups_TEMP<-df_groups_add%>%mutate(models.y=models.x)
        df_groups_TEMP<-df_groups_TEMP%>%mutate(RMSD=0)
        df_groups_TEMP<-unique(df_groups_TEMP)
        df_groups_add<-rbind(df_groups_add,df_groups_TEMP)
        
        df_groups_add<-df_groups_add%>%mutate(group=v_groups[j])
        df_groups_add<-df_groups_add%>%mutate(ligand_center=df_all$name[i])
        df_groups<-rbind(df_groups,df_groups_add)
      }
    }
    df_groups<-df_groups%>%group_by(ligand_center)%>%mutate(group_size=n())
    write.csv(df_groups,paste0("groups_fin/",df_all$name[i],".csv"),row.names = F)
  }
}

#copy pdb files to groups spb dir
if (!dir.exists("str")) {dir.create("str")}
if (!dir.exists("str_fin")) {dir.create("str_fin")}
i<-1
j<-1
k<-1
k<-2
#v_groups<-list.files(paste0("groups_fin/"))
for (j in 1:length(df_all$name)) {
  if(file.exists(paste0("groups_fin/",df_all$name[j],".csv"))){
    df_RMSD<-read.csv(paste0("groups_fin/",df_all$name[j],".csv"),stringsAsFactors = F)
    #  print(unique(df_RMSD$ligand_center))
    if(nrow(df_RMSD)>1){
      for (k in 1:nrow(df_RMSD)) {
        if (!dir.exists(paste0("str/",df_RMSD$ligand_center[k]))) { dir.create(paste0("str/",df_RMSD$ligand_center[k]))}
        if (!dir.exists(paste0("str/",df_RMSD$ligand_center[k],"/",df_RMSD$grop_number[k]))) {
          dir.create(paste0("str/",df_RMSD$ligand_center[k],"/",df_RMSD$grop_number[k]))}
        pdb<-read.pdb(paste0("pdb_second/",df_RMSD$ligand_center[k],"/",df_RMSD$models.y[k]))
        write.pdb(pdb,paste0("str/",df_RMSD$ligand_center[k],"/",df_RMSD$grop_number[k],"/",df_RMSD$models.y[k]))
      }
      df_RMSD<-df_RMSD%>%filter(models.y==models.x)
      for (k in 1:nrow(df_RMSD)) {
        pdb<-read.pdb(paste0("pdb_second/",df_RMSD$ligand_center[k],"/",df_RMSD$models.y[k]))
        write.pdb(pdb,paste0("str_fin/",df_RMSD$ligand_center[k],"_",df_RMSD$grop_number[k],"_",df_RMSD$models.y[k]))
      }
    }
  }
}
#energy bonding
df_all<-read.csv(paste0(part_name,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
i<-1
if (!dir.exists("log_fin")) {dir.create("log_fin")}
#if (!dir.exists("plot")) {dir.create("plot")}
df_log<-read.csv("df_log_all.csv",stringsAsFactors = F)
df_log<-df_log%>%mutate(models.y=paste0("frame_",new_number,".pdb"))
for (i in 1:nrow(df_all)) {
  if(file.exists(paste0("groups_fin/",df_all$name[i],".csv"))){
    df_groups<-read.csv(paste0("groups_fin/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_fin<-left_join(df_groups,df_log,by=c("models.y","ligand_center"="name"))
    write.csv(df_fin,paste0("log_fin/",df_all$name[i],".csv"),row.names = F)
  }
}

df_fin<-read.csv(paste0("log_fin/",df_all$name[1],".csv"),stringsAsFactors = F)
for (i in 2:length(df_all$name)) {
  if(file.exists(paste0("groups_fin/",df_all$name[i],".csv"))){
    df_fin_add<-read.csv(paste0("log_fin/",df_all$name[i],".csv"),stringsAsFactors = F)
    df_fin<-rbind(df_fin,df_fin_add)
  }
}

p<-ggplot(df_fin)+
  geom_freqpoly(aes(x=affinity,colour=group),binwidth=0.3)+
  facet_grid(receptor~ligand)+
  theme_bw()+guides(color = "none", size = "none")
ggsave(p,filename = paste0("ligand_energy.png"), width = 30, height = 20, units = c("cm"), dpi = 200 ) 
write.csv(df_fin,"log_fin.csv",row.names = F)