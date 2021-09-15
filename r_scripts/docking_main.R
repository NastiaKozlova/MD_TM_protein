part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)

v_ligand<-list.files(paste0("start/ligand_start/"))
a<-c()
for (i in 1:length(v_ligand)) {
    b<-strsplit(v_ligand[i],split = ".",fixed = T)[[1]][1]
    a<-c(a,b)
}
v_ligand<-a
  
df_ligand<-data.frame(matrix(ncol=2,nrow=length(v_ligand)))
colnames(df_ligand)<-c("ligand","c")
df_ligand$ligand<-v_ligand
  
df_active_center<-read.csv(paste0("start/active_center.csv"),stringsAsFactors = F)
v_center<-unique(df_active_center$type)
df_center<-data.frame(matrix(ncol=2,nrow=length(v_center)))
colnames(df_center)<-c("center","c")
df_center$center<-v_center
  
df_ligand_center<-full_join(df_ligand,df_center,by="c")
df_ligand_center$c<-NULL
write.csv(df_ligand_center,paste0("start/ligand_center.csv"),row.names = F)

#docking
if (!dir.exists(paste0(part_start,"MD_analysis/docking/"))){dir.create(paste0(part_start,"MD_analysis/docking/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/docking/docking_first/"))){dir.create(paste0(part_start,"MD_analysis/docking/docking_first/"))}
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")
#if(!dir.exists(paste0(part,name,"/docking"))){dir.create(paste0(part,name,"/docking"))}
if(!dir.exists(part_name)){dir.create(part_name)}
system(command = paste0("cp -r ",part_start,"/start/ligand_start/ ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("cp ",part_start,"/start/active_center.csv ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("cp ",part_start,"/start/ligand_center.csv ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"MD_analysis/docking/receptor_start/ ",part_name),ignore.stdout=T,wait = T)
df_ligand_center<-read.csv(paste0(part_name,"ligand_center.csv"),stringsAsFactors = F)
df_ligand_center<-df_ligand_center%>%mutate(c="C")
v_receptor<-list.files(paste0(part_name,"receptor_start/"))
a<-c()
for (i in 1:length(v_receptor)) {
  b<-strsplit(v_receptor[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_receptor<-a
  
df_receptor<-data.frame(matrix(ncol=2,nrow=length(v_receptor)))
colnames(df_receptor)<-c("receptor","c")
df_receptor$receptor<-v_receptor
df_receptor<-df_receptor%>%mutate(c="C")
df_all<-full_join(df_receptor,df_ligand_center,by="c")
df_all$c<-NULL
write.csv(df_all,paste0(part_name,"df_all.csv"),row.names = F)
if (!dir.exists(paste0(part_name,"ligand/"))){dir.create(paste0(part_name,"ligand/"))}
if (!dir.exists(paste0(part_name,"receptor/"))){dir.create(paste0(part_name,"receptor/"))}
i<-1
for (i in 1:nrow(df_receptor)) {
  system(command = paste0(part_start,"programs/MGLTools-1.5.7/bin/pythonsh ",part_start,"programs/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ",
                          part_name,"receptor_start/",df_receptor$receptor[i],".pdb -o ",part_name,"receptor/",df_receptor$receptor[i],".pdbqt ",
                          "-A None"))
}

for (i in 1:nrow(df_ligand)) {
  system(command = paste0("obabel ",part_name,"ligand_start/",df_ligand$ligand[i], ".pdb -O ",part_name,"ligand/",df_ligand$ligand[i], ".pdbqt"),ignore.stdout=T,wait = T)
}


system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_script.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_name,"script_fin.txt "),ignore.stdout=T,wait = T)
system(command = paste0(part_name,"script_fin.txt"),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/convert_pdbqt_to_pdb.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_name,"convert_pdbqt_to_pdb.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_name,"convert_pdbqt_to_pdb.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_pre_analysis.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_group_structure.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_interactions.R ",part_name),ignore.stdout=T,wait = T)
#part_start<-paste0(part_start,"MD_analysis/")
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/analysis.R ",part_name),ignore.stdout=T,wait = T)

df_fin<-read.csv(paste0(part_name,"din/log_fin.csv"),stringsAsFactors =  F)
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors =  F)
df_all_systems<-df_all_systems%>%mutate(system=paste(Structure,Membrane))
df_all_systems<-df_all_systems%>%mutate(system_name=paste0("charmm-gui-",system_name))
df_fin<-left_join(df_fin,df_all_systems,by=c("receptor"="system_name"))
p<-ggplot(df_fin)+
  geom_freqpoly(aes(x=affinity,colour=group),binwidth=0.3)+
  facet_grid(ligand~system)+
  theme_bw()+guides(colour = "none")
ggsave(p,filename = paste0(part_start,"plot/ligand_energy.png"), width = 20, height = 15, units = c("cm"), dpi = 200 ) 
p<-ggplot(df_fin)+
  geom_freqpoly(aes(x=affinity,colour=group),binwidth=0.3)+
  facet_grid(system~ligand)+
  theme_bw()+guides(colour = "none")
ggsave(p,filename = paste0(part_start,"plot/ligand_energy_reverce.png"), width = 20, height = 15, units = c("cm"), dpi = 200 ) 
