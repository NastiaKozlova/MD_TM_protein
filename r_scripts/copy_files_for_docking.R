part_name <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
part_TEMP<-strsplit(part_name,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_search<-part_TEMP[2]

setwd(part_start)
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")
if(!dir.exists(part_name)){dir.create(part_name)}
part_name<-paste0(part_name,v_search,"/")
#if(!dir.exists(part_name)){dir.create(part_name)}
if(!dir.exists(part_name)){dir.create(part_name)}
print("test")
#if(!dir.exists(part_name)){dir.create(part_name)}
if(!dir.exists(paste0(part_name,"ligand"))){dir.create(paste0(part_name,"ligand"))}
system(command = paste0("cp -r ",part_start,"/start/ligand_start/ ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("cp -r ",part_start,"programs ",part_name),ignore.stdout=T,wait = T)
#ligands_prepare
v_ligand<-list.files(paste0(part_start,"start/ligand_start/"))
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
if(v_search=="center"){
    system(command = paste0("cp ",part_start,"start/active_center.csv ", part_name),ignore.stdout=T,wait = T)
}
