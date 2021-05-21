part_start <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_start)
part<-part_start

#library(dplyr)
#library(bio3d)
#library(ggplot2)


v_first_bond<-list.files(paste0(part_start,"docking/docking_first/din/interaction_fin/"))
df_first_bond_start<-read.csv(paste0(part_start,"docking/docking_first/din/interaction_fin/",v_first_bond[1]),stringsAsFactors = F)
for (i in 2:length(v_first_bond)) {
  df_first_bond_add<-read.csv(paste0(part_start,"docking/docking_first/din/interaction_fin/",v_first_bond[i]),stringsAsFactors = F)
  df_first_bond_start<-rbind(df_first_bond_start,df_first_bond_add)
  df_first_bond_start<-df_first_bond_start%>%filter(number_interactions>80)
}
df_first_bond_start<-df_first_bond_start%>%filter(number_interactions>80)

df_first_bond_start<-df_first_bond_start%>%select(resno,number_interactions, system, center, ligand, grops, grops_number, 
                                                  persent_interactions)
df_first_bond_start<-df_first_bond_start%>%select(resno, system, center, ligand, grops, grops_number, 
                                                  persent_interactions)
df_first_bond_start<-df_first_bond_start%>%mutate(name=paste0(system,"_", ligand,"_",center ))

v_second_bond<-list.files(paste0(part_start,"docking/docking_second/din/interaction_fin/"))
df_second_bond_start<-read.csv(paste0(part_start,"docking/docking_second/din/interaction_fin/",v_second_bond[1]),stringsAsFactors = F)
for (i in 2:length(v_second_bond)) {
  df_second_bond_add<-read.csv(paste0(part_start,"docking/docking_second/din/interaction_fin/",v_second_bond[i]),stringsAsFactors = F)
  df_second_bond_start<-rbind(df_second_bond_start,df_second_bond_add)
  df_second_bond_start<-df_second_bond_start%>%filter(number_interactions>80)
}
df_second_bond_start<-df_second_bond_start%>%select(resno,number_interactions, system,
                                                    center, ligand, grops, grops_number, persent_interactions)
df_second_bond_start<-df_second_bond_start%>%filter(number_interactions>80)
df_second_bond_start<-df_second_bond_start%>%select(resno,system, center, ligand, grops, grops_number, persent_interactions)
df_fin_semi<-semi_join(df_first_bond_start,df_second_bond_start,by=c("center", "ligand","resno"))
df_fin_anti<-anti_join(df_first_bond_start,df_second_bond_start,by=c("center", "ligand","resno"))
df_fin<-left_join(df_first_bond_start,df_second_bond_start,by=c("center", "ligand","name"="system"))
write.csv(df_fin,paste0(part_start,"docking/df_interactions.csv"),row.names = F)



df_fin_no<-df_fin%>%filter(resno.x!=resno.y)
df_fin_yes<-df_fin%>%filter(resno.x==resno.y)
df_first_bond_start$name<-NULL
df_first_bond_start<-df_first_bond_start%>%mutate(ligand_center=paste0(ligand,"_",center))
v_ligand_center<-unique(df_first_bond_start$ligand_center)
df_second_bond_start<-df_second_bond_start%>%mutate(ligand_center=paste0(ligand,"_",center))
df_second_bond_start<-df_second_bond_start%>%mutate(tes="second")
i<-1
for (i in 1:length(v_ligand_center)) {
  df_first_bond_TEMP<-system%>%filter(ligand_center==v_ligand_center[i])
  df_second_bond_start$tes[df_second_bond_start$resno%in%df_first_bond_TEMP$resno&
                           df_second_bond_start$ligand_center==v_ligand_center[i]]<-"both"
}

df_second_bond_no<-df_second_bond_start%>%filter(tes=="both")
df_second_bond_no<-df_second_bond_no%>%mutate(system="start")
df_second_bond_no<-unique(df_second_bond_no)
df_second_bond_no<-df_second_bond_no%>%mutate(temp=paste0(resno,system,center,ligand,ligand_center))
df_first_bond_start<-df_first_bond_start%>%mutate(tes="first")
df_first_bond_start<-df_first_bond_start%>%mutate(temp=paste0(resno,system,center,ligand,ligand_center))
#df_first_bond_start<-df_first_bond_start%>%mutate(tes="first")

df_first_bond_start$tes[df_first_bond_start$temp%in%df_second_bond_no$temp]<-NA
df_first_bond_start<-df_first_bond_start%>%filter(!is.na(tes))
df_first_bond_start$temp<-NULL
df_fin<-rbind(df_first_bond_start,df_second_bond_start)

df_fina<-df_fin%>%filter(tes=="first")
df_finb<-df_fin%>%filter(tes=="both")
df_finc<-df_fin%>%filter(tes=="second")



part_start<-paste0(part_start,"docking/docking_second/")
setwd(part_start)
df_groups_start<-read.csv("din/df_groups_start.csv",stringsAsFactors = F)
df_all<-read.csv("df_all.csv",stringsAsFactors = F)
df_all<-df_all%>%mutate(name=paste0(receptor,"_",ligand,"_",center))
df_groups_start<-full_join(df_groups_start,df_all,by=c("ligand_center"="name"))
