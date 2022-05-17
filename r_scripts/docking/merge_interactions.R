part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_name)
setwd("din")


df_all<-read.csv(paste0("df_merge_structure_log.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center))

#if (dir.exists(paste0("interaction/"))) { system(command = paste0("rm -r ",part_name,"din/interaction/"))}
#if (dir.exists(paste0("interaction_TEMP/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_TEMP/"))}
if (dir.exists(paste0("interaction_serf/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_serf/"))}

if (!dir.exists(paste0("interaction_serf/"))) { dir.create(paste0("interaction_serf/"))}
i<-1
j<-1
p<-1
v_structure<-unique(df_all$name.x)
for (j in 1:length(v_structure)) {
  df_complex<-df_all%>%filter(name.x==v_structure[i])
  pdb<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[j],".pdb"))
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%mutate(number_interactions=0)
  df_pdb<-df_pdb%>%mutate(tested_structure=0)
  df_pdb<-df_pdb%>%filter(total_structure=nrow(df_complex))
  for (p in 1:nrow(df_complex)) {
    if(file.exists(paste0("interaction/",df_all$receptor_ligand[p],"/",df_complex$new_number[p],".csv"))){
      df_protein<-read.csv(paste0("interaction/",df_all$receptor_ligand[p],"/",df_complex$new_number[p],".csv"),
                           stringsAsFactors = F) 
      df_pdb$number_interactions[df_pdb$resno%in%df_protein$resno]<-df_pdb$number_interactions[df_pdb$resno%in%df_protein$resno]+1
      df_pdb$tested_structure<-df_pdb$total_structure+1
      
    }
  }
  write.csv(df_protein,
            paste0("interaction_serf/",df_all$receptor_ligand[j],"/",df_all$new_number[j],".csv"),
            row.names = F)
}


#i<-1
#j<-1
#df_topology<-df_all%>%select(receptor_ligand,receptor,ligand, center,size_of_group)
#df_topology<-unique(df_topology)
#for (j in 1:nrow(df_topology)) {
#  if (!dir.exists(paste0("interaction_TEMP/",df_topology$receptor_ligand[j]))){dir.create(paste0("interaction_TEMP/",df_topology$receptor_ligand[j]))}
#  pdb<-read.pdb(paste0(part_name,"receptor_start/",df_topology$receptor[j],".pdb"))
#  df_pdb<-pdb$atom
#  df_pdb<-df_pdb%>%filter(elety=="CA")
#  df_pdb<-df_pdb%>%select(type,resid,resno, x,y,z)
#  df_pdb<-df_pdb%>%mutate(number_interactions=0)
#  df_pdb<-df_pdb%>%mutate(receptor_ligand=df_topology$receptor_ligand[j])
#  df_pdb<-df_pdb%>%mutate(receptor=df_topology$receptor[j])
#  df_pdb<-df_pdb%>%mutate(center=df_topology$center[j])
#  df_pdb<-df_pdb%>%mutate(ligand=df_topology$ligand[j])
#  df_pdb<-df_pdb%>%mutate(size_of_group=df_topology$size_of_group[j])
#  v_frame<-list.files(paste0("interaction/",df_topology$receptor_ligand[j],"/",df_topology$size_of_group[j],"/",df_topology$center[j],"/"))
#  for (q in 1:length(v_frame)) {
#    df_interaction<-read.csv(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j],"/",df_all$center[j],"/",df_all$models.y[j],".csv"),stringsAsFactors = F)
#    colnames(df_interaction)<-c(colnames(df_interaction)[2],colnames(df_interaction)[1])
#    df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]<-df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]+1
#  }
#  df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/length(v_frame)*100)
#  write.csv(df_pdb,paste0("interaction_TEMP/",df_topology$receptor_ligand[j],"/",df_topology$size_of_group[j],"_",df_topology$center[j],".csv"),row.names = F)
#}
#i<-1
#j<-1#

#v_system<-list.files(paste0("interaction_TEMP/"))
#for (i in 1:length(v_system)) {
#  v_frame<-list.files(paste0("interaction_TEMP/",v_system[i]))
#  df_pdb<-read.csv(paste0("interaction_TEMP/",v_system[i],"/",v_frame[1]))
#  for (q in 2:length(v_frame)) {
#    df_pdb_add<-read.csv(paste0("interaction_TEMP/",v_system[i],"/",v_frame[q]))
#    df_pdb<-rbind(df_pdb,df_pdb_add)
#  }
#  df_pdb<-df_pdb%>%mutate(sort=paste(resid,resno,receptor_ligand,size_of_group))
#  df_pdb<-df_pdb%>%group_by(sort)%>%mutate(total_interactions=sum(number_interactions))
#  df_pdb<-ungroup(df_pdb)
#  df_pdb<-df_pdb%>%mutate(total_persent_interactions=total_interactions/size_of_group*100)
#  df_pdb<-df_pdb%>%select(resid,resno,
#                          receptor_ligand,receptor,ligand,size_of_group,             
#                          total_interactions,total_persent_interactions)
#  df_pdb<-unique(df_pdb)
#  write.csv(df_pdb,paste0("interaction_fin/",v_system[i],".csv"),row.names = F)
#}
#v_groups<-list.files(paste0("interaction_fin/"))
#df_pdb<-read.csv(paste0("interaction_fin/",v_groups[1]),stringsAsFactors =  F)
#if(length(v_groups)>1){
#  for (i in 2:length(v_groups)) {
#    df_pdb_add<-read.csv(paste0("interaction_fin/",v_groups[i]),stringsAsFactors =  F)
#    df_pdb<-rbind(df_pdb,df_pdb_add)
#  }
#}

#df_pdb<-df_pdb%>%filter(total_persent_interactions>0)

#p<-ggplot(data=df_pdb)+geom_freqpoly(aes(x=total_persent_interactions))+theme_bw()+facet_grid(size_of_group~receptor_ligand)
#ggsave(p,filename = paste0(part_name,"interaction_ligand_receptor.png"), width = 12, height = 12, units = c("cm"), dpi = 1000 ) 
