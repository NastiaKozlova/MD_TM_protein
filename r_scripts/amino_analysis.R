part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(dplyr)
library(bio3d)
test_10<-seq(from=0,to=1000,by=10)
df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors =  F)
part<-paste0(part_start,"MD_analysis/")
setwd(part)
parta<-paste0(part_start,"MD/")
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)

df_all_systems$system_name<-as.character(df_all_systems$system_name)

df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))

for (q in 1:nrow(df_all_systems)) {
  if(file.exists(paste0("fin_data/str_data/",df_all_systems$fin_name[q],".csv"))){
    if(file.exists(paste0("fin_data/aminoacids_importance/",df_all_systems$fin_name[q],".csv"))){
      df_start<-read.csv(paste0("fin_data/str_data/",df_all_systems$fin_name[q],".csv"),stringsAsFactors = F)
      df_add<-read.csv(paste0("fin_data/aminoacids_importance/",df_all_systems$fin_name[q],".csv"),stringsAsFactors = F)
#      df_start<-df_start%>%filter(ligand=="D-Lactose")
#      df_start<-df_start%>%filter(grops_number==1)
      df_start<-left_join(df_start,df_add,by="resno")
#      df_start<-df_start[(df_start$persent_interactions>80|df_start$importance>1),]
      df_start<-df_start%>%mutate(system=df_all_systems$fin_name[q])
      df_start<-left_join(df_start,df_all_systems,by=c("system"="fin_name"))
      write.csv(df_start,paste0("fin_data/seq_analysis/",df_all_systems$fin_name[q],".csv"),row.names = F)
    }
  }
}

df_start<-df_start%>%filter(resno<0)
for (q in 1:nrow(df_all_systems)) {
  if(file.exists(paste0("fin_data/seq_analysis/",df_all_systems$fin_name[q],".csv"))){
    df_start_add<-read.csv(paste0("fin_data/seq_analysis/",df_all_systems$fin_name[q],".csv"),stringsAsFactors = F)
    df_start<-rbind(df_start,df_start_add)
  }
}
write.csv(df_start,paste0("fin_data/seq_analysis.csv"),row.names = F)   
