part_start <- commandArgs(trailingOnly=TRUE)
if(!dir.exists(paste0(part_start,"MD_analysis/docking/docking_first/"))){dir.create(paste0(part_start,"MD_analysis/docking/docking_first/"))}

library(dplyr)
library(bio3d)
#library(readr)
filter_structure<-function(df_pdb,x_min,x_max,y_min,y_max,z_min,z_max){
  df_pdb_filtered<-df_pdb%>%filter(x>x_min)
  df_pdb_filtered<-df_pdb_filtered%>%filter(x<x_max)
  df_pdb_filtered<-df_pdb_filtered%>%filter(y>y_min)
  df_pdb_filtered<-df_pdb_filtered%>%filter(y<y_max)
  df_pdb_filtered<-df_pdb_filtered%>%filter(z>z_min)
  df_pdb_filtered<-df_pdb_filtered%>%filter(z<z_max)
  if(nrow(df_pdb_filtered)>3){
    df_pdb_filtered<-df_pdb_filtered%>%select(type,resid,resno)
    colnames(df_pdb_filtered)<-c("type","amino","resno") 
    df_pdb_filtered<-df_pdb_filtered%>%mutate(type=df_structure$type[i])
  }
  return(df_pdb_filtered)
}
setwd(part_start)
df_all_systems<-read.csv("start/all_systems.csv",stringsAsFactors = F)
part_name<-paste0(part_start,"MD_analysis/")
setwd(part_name)

j<-1
if(!dir.exists(paste0("docking/active_center_TEMP/"))){dir.create("docking/active_center_TEMP/")}
for (j in 1:nrow(df_all_systems)) {
  start<-read.pdb(paste0("docking/receptor_start/charmm-gui-",df_all_systems$system_name[j],".pdb"))
  df_hbonds<-read.csv(paste0("din/charmm-gui-",df_all_systems$system_name[j],"/hbonds_8.csv"),stringsAsFactors = F)
  
  df_pdb<-start$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-left_join(df_pdb,df_hbonds,by=c("resno"="number","resid"="amino" )) 
  df_pdb<-df_pdb%>%filter(persent>50)
  x_min<-round(min(df_pdb$x),digits = 0)
  y_min<-round(min(df_pdb$y),digits = 0)
  z_min<-round(min(df_pdb$z),digits = 0)
  x_max<-round(max(df_pdb$x),digits = 0)
  y_max<-round(max(df_pdb$y),digits = 0)
  z_max<-round(max(df_pdb$z),digits = 0)
  df_structure<-data.frame(matrix(ncol=7,nrow = 0))
  colnames(df_structure)<-c("type","x_min","y_min","z_min","x_max","y_max","z_max")
  x_len<-seq(from=x_min,to=x_max,by=10)
  y_len<-seq(from=y_min,to=y_max,by=10)
  z_len<-seq(from=z_min,to=z_max,by=10)
  x<-x_len[1]
  y<-y_len[1]
  z<-z_len[1]
  for (x in x_len) {
    for (y in y_len) {
      for (z in z_len) {
        df_structure_add<-data.frame(matrix(ncol=7,nrow = 1))
        colnames(df_structure_add)<-c("type","x_min","y_min","z_min","x_max","y_max","z_max")
        df_structure_add$type<-paste0("serf_x_",x,"_y_",y,"_z_",z)
        df_structure_add$x_min<-x-10
        df_structure_add$x_max<-x+10
        df_structure_add$y_min<-y-10
        df_structure_add$y_max<-y+10
        df_structure_add$z_min<-z-10
        df_structure_add$z_max<-z+10
        df_structure<-rbind(df_structure,df_structure_add)
      }
    }
  }
  df_structure<-unique(df_structure)
  i<-1
  df_type<-data.frame(matrix(ncol=3,nrow = 0))
  colnames(df_type)<-c("type","amino","resno") 
  for (i in 1:nrow(df_structure)) {
    df_pdb_filtered<-filter_structure(df_pdb = df_pdb,
                                      x_min = df_structure$x_min[i],
                                      x_max = df_structure$x_max[i],
                                      y_min = df_structure$x_min[i],
                                      y_max = df_structure$x_max[i],
                                      z_min = df_structure$x_min[i],
                                      z_max = df_structure$x_max[i])
    if(ncol(df_pdb_filtered)==3){
      df_type<-rbind(df_type,df_pdb_filtered)
      
    }else{df_structure$type[i]<-NA}
  }
  df_structure<-df_structure%>%filter(!is.na(type))
  df_type<-read.csv(paste0(part_start,"start/active_center.csv"),stringsAsFactors = F)
  df_type<-unique(df_type)
  df_type<-df_type%>%mutate(receptor=paste0("charmm-gui-",df_all_systems$system_name[j]))
  write.csv(df_type,paste0("docking/active_center_TEMP/charmm-gui-",df_all_systems$system_name[j],"_active_center.csv"),row.names = F)
}
j<-3
df_type<-read.csv(paste0("docking/active_center_TEMP/charmm-gui-",df_all_systems$system_name[1],"_active_center.csv"),stringsAsFactors =  F)
print(nrow(df_type))
for (j in 2:nrow(df_all_systems)) {
  df_type_add<-read.csv(paste0("docking/active_center_TEMP/charmm-gui-",df_all_systems$system_name[j],"_active_center.csv"),stringsAsFactors =  F)
  print(paste(df_all_systems$system_name[j],nrow(df_type_add)))
  df_type<-rbind(df_type,df_type_add)
  print(nrow(df_type))
}

write.csv(df_type,paste0("docking/docking_first/active_center.csv"),row.names = F)
