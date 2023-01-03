part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)

part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")

if(!dir.exists(part_name)){dir.create(part_name)}
if(!dir.exists(paste0(part_name,"ligand"))){dir.create(paste0(part_name,"ligand"))}
system(command = paste0("cp -r ",part_start,"/start/ligand_start/ ",part_name),ignore.stdout=T,wait = T)
#ligands_prepare
v_ligand<-list.files(paste0("start/ligand_start/"))
a<-c()
for (i in 1:length(v_ligand)) {
  b<-strsplit(v_ligand[i],split = ".",fixed = T)[[1]][1]
  a<-c(a,b)
}
v_ligand<-a
df_ligand<-data.frame(matrix(ncol=2,nrow=length(v_ligand)))
colnames(df_ligand)<-c("ligand","C")
df_ligand$ligand<-v_ligand
for (i in 1:nrow(df_ligand)) {
  system(command = paste0("obabel ",part_name,"ligand_start/",df_ligand$ligand[i], ".pdb -O ",part_name,"ligand/",df_ligand$ligand[i], ".pdbqt"),ignore.stdout=T,wait = T)
}

#prepare_receptor

if (!dir.exists(paste0(part_name,"receptor/"))){dir.create(paste0(part_name,"receptor/"))}
if (!dir.exists(paste0(part_name,"active_center/"))){dir.create(paste0(part_name,"active_center/"))}
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)
system(command = paste0("cp -r ",part_start,"/MD_analysis/docking/receptor_start/ ",part_name),ignore.stdout=T,wait = T)
df_receptor<-df_all_systems%>%mutate(receptor=paste0("charmm-gui-",system_name))
df_receptor<-df_receptor%>%mutate(c="C")

for (i in 1:nrow(df_receptor)) {
  system(command = paste0(part_start,"programs/MGLTools-1.5.7/bin/pythonsh ",part_start,"programs/MGLTools-1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r ",
                          part_name,"receptor_start/",df_receptor$receptor[i],".pdb -o ",part_name,"receptor/",df_receptor$receptor[i],".pdbqt ",
                          "-A None"))
}
j<-1
df_active_center<-read.csv(paste0(part_name,"active_center.csv"),stringsAsFactors = F)
df_active_center<-df_active_center%>%mutate(C=NA)

df_ligand_center<-left_join(df_active_center,df_ligand,by="C")
df_ligand_center<-df_ligand_center%>%select(type,receptor,ligand)
df_ligand_center<-unique(df_ligand_center)
write.csv(df_ligand_center,paste0(part_name,"ligand_center.csv"),row.names =  F)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_script.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/convert_pdbqt_to_pdb.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_name,"convert_pdbqt_to_pdb.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_name,"convert_pdbqt_to_pdb.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/convert_log_to_csv.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("chmod +x ",part_name,"prepare_log_csv.py "),ignore.stdout=T,wait = T)
system(command = paste0("python3 ", part_name,"prepare_log_csv.py"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_pre_analysis.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_interactions.R ",part_name),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/RMSD_group_structure.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/calibration_docking_group_structure.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"docking_group_structure.R ",part_analysis,",",1),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"RMSD_merge_docking_parts.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"calibration_merge_structure.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_docking_parts.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"merge_interactions.R ",part_analysis),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_scriprs,"complex_structure_surf.R ",part_analysis),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_scriprs,"atom_interactions_surf.R ",part_analysis),ignore.stdout=T,wait = T)
