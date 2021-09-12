part_start <- commandArgs(trailingOnly=TRUE)
#part_start<-part_name
setwd(part_start)
library(bio3d)
library(readr)
library(dplyr)
library(ggplot2)
v_rmsd<-4


#library(dplyr)
#library(bio3d)
#library(ggplot2)

#df_RMSD<-read.csv(paste0("docking/df_RMSD.csv"),stringsAsFactors = F)
v_first_bond<-list.files(paste0("din/interaction_fin/"))
df_first_bond_start<-read.csv(paste0("din/interaction_fin/",v_first_bond[1]),stringsAsFactors = F)
for (i in 2:length(v_first_bond)) {
  df_first_bond_add<-read.csv(paste0("din/interaction_fin/",v_first_bond[i]),stringsAsFactors = F)
  df_first_bond_start<-rbind(df_first_bond_start,df_first_bond_add)
  df_first_bond_start<-df_first_bond_start%>%filter(persent_interactions>80)
}
df_first_bond_start<-df_first_bond_start%>%filter(persent_interactions>80)
df_first_bond_start<-df_first_bond_start%>%filter(persent_interactions==100)

#df_first_bond_start<-df_first_bond_start[df_first_bond_start$system%in%df_RMSD$fin_model,]
df_first_bond_start<-df_first_bond_start%>%mutate(receptor=system)
#for (i in 1:nrow(df_first_bond_start)) {
#  df_first_bond_start$receptor[i]<-strsplit(df_first_bond_start$system[i],split = "_",fixed = T)[[1]][1]  
#}
#df_first_bond_start<-df_first_bond_start%>%select(receptor,ligand, center, resid,resno, x,y,z,system,grops,grops_number)
df_first_bond_start<-df_first_bond_start%>%select(receptor,ligand, center, resid,resno, grops,grops_number,system)
df_first_bond<-unique(df_first_bond_start)


df_first_bond<-df_first_bond%>%mutate(test_complex_amino=paste(receptor,ligand, center, resid,resno,sep="_"))
df_first_bond<-df_first_bond%>%group_by(test_complex_amino)%>%mutate(number_interactons_complex=n())
df_first_bond<-ungroup(df_first_bond)

df_first_bond_system<-df_first_bond%>%mutate(test_complex=paste(receptor,ligand, center, sep="_"))
df_first_bond_system<-df_first_bond_system%>%select(receptor,ligand,center,test_complex,system,grops)
df_first_bond_system<-unique(df_first_bond_system)

df_first_bond_system<-df_first_bond_system%>%group_by(test_complex)%>%mutate(number_models_complex=n())
df_first_bond_system<-ungroup(df_first_bond_system)
df_first_bond_system<-df_first_bond_system%>%select(receptor, ligand,center,number_models_complex)
#df_first_bond_system<-df_first_bond_system%>%group_by(test_system)%>%mutate(number_models_system=n())
df_first_bond_system<-unique(df_first_bond_system)

df_first_bond<-ungroup(df_first_bond)
df_first_bond_system<-ungroup(df_first_bond_system)
df_first_bond_system$system<-NULL
df_first_bond$system<-NULL
#df_first_bond_system<-ungroup(df_first_bond_system)
#df_first_bond_system<-unique(df_first_bond_system)


df_first_bonda<-left_join(df_first_bond,df_first_bond_system,
                          by = c("receptor", "ligand", "center"))

df_first_bonda<-df_first_bonda%>%mutate(occurence_interctions_complex=number_interactons_complex/number_models_complex*100)

df_first_bonda<-df_first_bonda%>%select(receptor, ligand, center, resid, resno, occurence_interctions_complex)
df_first_bonda<-unique(df_first_bonda)
df_first_bonda<-df_first_bonda%>%filter(occurence_interctions_complex>50)
p<-ggplot(data = df_first_bonda)+
  geom_text(aes(x=resno,y=resid,colour=occurence_interctions_complex,label=resno))+facet_grid(receptor~ligand)+
  theme_bw()+guides(colour="none")
ggsave(p,   filename = paste0("docking_aminoacid_interaction.png"), width = 40, height = 20, units = c("cm"), dpi = 200 ) 


write.csv(df_first_bonda,"docking_aminoacid_interactions.csv",row.names = F)
df_topology<-read.csv("din/df_topology.csv",stringsAsFactors = F)
df_fin<-read.csv("din/log_fin.csv",stringsAsFactors = F)
df_fin<-df_fin%>%mutate(receptor_fin=NA)
for (i in 1:nrow(df_fin)) {
  df_fin$receptor_fin[i]<-strsplit(df_fin$receptor[i],split = "_",fixed = T)[[1]][1]
}
p<-ggplot(df_fin)+
  geom_freqpoly(aes(x=affinity,colour=grop_number),binwidth=0.3)+
  facet_grid(receptor_fin~ligand)+
  theme_bw()
ggsave(p,filename = paste0("ligand_energy.png"), width = 30, height = 20, units = c("cm"), dpi = 200 ) 
