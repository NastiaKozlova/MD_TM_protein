part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-4

setwd(part_name)

main_part<-list.files("interaction")

df_all<-read.csv(paste0(part_name,"din/df_merge_structure_log.csv"),stringsAsFactors = F)


j<-1
setwd("din")
if (!dir.exists(paste0("interaction/"))) { dir.create(paste0("interaction/"))}
i<-1
for (j in 1:nrow(df_all)) {
  if(file.exists(paste0("groups_fin/",df_all$ligand_center[j],".csv"))){
    
    a<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[i],".pdb"))
    b<-read.pdb(paste0("pdb_second/",df_all$ligand_center[j],"/",df_all$models.y[j]))
    bs<-binding.site(a,b)
    m<-bs$resnames
    a<-c()
    b<-c()
    y<-1
    for (y in 1:length(m)) {
      p<-strsplit(m[y],split = " ",fixed = T)[[1]][2]
      a<-c(a,p)
      p<-strsplit(m[y],split = " ",fixed = T)[[1]][1]
      b<-c(b,p)
    }
    a<-as.numeric(a)
    df_protein<-data.frame(matrix(ncol=2,nrow=length(a)))
    colnames(df_protein)<-c("resid","resno")
    df_protein$resid<-a
    df_protein$resno<-b
    if (!dir.exists(paste0("interaction/",df_all$receptor_ligand[j]))) { dir.create(paste0("interaction/",df_all$receptor_ligand[j]))}
    if (!dir.exists(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j]))) {
      dir.create(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j]))}
    if (!dir.exists(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j],"/",df_all$size_of_group[j],"/"))) {
      dir.create(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j],"/",df_all$center[j],"/"))}
    write.csv(df_protein,
              paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j],"/",df_all$center[j],"/",df_all$models.y[j],".csv"),
              row.names = F)
  }
}

if (!dir.exists("din/interaction_fin/")){dir.create("din/interaction_fin/")}
i<-1
j<-1
df_topology<-read.csv("din/df_topology.csv",stringsAsFactors = F)
v_name<-list.files(paste0("din/interaction/"))
df_topology<-df_topology[df_topology$name%in%v_name,]
df_topology<-df_topology
df_topology$name_log<-NULL
df_topology$run<-NULL
df_topology<-unique(df_topology)
for (j in 1:nrow(df_topology)) {
  pdb<-read.pdb(paste0("receptor_start/",df_topology$receptor[j],".pdb"))
  v_groups<-list.files(paste0("din/interaction/",df_topology$name[j]))
  if(length(v_groups)>0){
    print(j)
    for (i in 1:length(v_groups)) {
      v_frame<-list.files(paste0("din/interaction/",df_all$receptor_ligand[j]))
      df_pdb<-pdb$atom
      df_pdb<-df_pdb%>%filter(elety=="CA")
      df_pdb<-df_pdb%>%mutate(number_interactions=0)
      df_pdb<-df_pdb%>%mutate(system=df_topology$receptor[j])
      df_pdb<-df_pdb%>%mutate(center=df_topology$center[j])
      df_pdb<-df_pdb%>%mutate(ligand=df_topology$ligand[j])
      df_pdb<-df_pdb%>%mutate(grops=v_groups[i])
      df_pdb<-df_pdb%>%mutate(grops_number=i)
      for (q in 1:length(v_frame)) {
        df_interaction<-read.csv(paste0("interaction/",df_all$receptor_ligand[j],"/",df_all$size_of_group[j],"_",df_all$center[j],"/",df_all$models.y[j],".csv"),stringsAsFactors = F)
        colnames(df_interaction)<-c(colnames(df_interaction)[2],colnames(df_interaction)[1])
        df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]<-df_pdb$number_interactions[df_pdb$resno%in%df_interaction$resno]+1
      }
      df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/length(v_frame)*100)
      write.csv(df_pdb,paste0("din/interaction_fin/",df_topology$name[j],"_",v_groups[i],".csv"),row.names = F)
    }
  } 
}
v_groups<-list.files(paste0("din/interaction_fin/"))
df_pdb<-read.csv(paste0("din/interaction_fin/",v_groups[1]),stringsAsFactors =  F)
for (i in 2:length(v_groups)) {
  df_pdb_add<-read.csv(paste0("din/interaction_fin/",v_groups[i]),stringsAsFactors =  F)
  df_pdb<-rbind(df_pdb,df_pdb_add)
}

df_pdb<-df_pdb%>%mutate(aminoacids=paste(resno,resid))
df_pdb<-df_pdb%>%select(resno, x, y, z, number_interactions, system, center, ligand, 
                        grops, grops_number,persent_interactions, aminoacids )
write.csv(df_pdb,"interaction_fin.csv",row.names = F)
