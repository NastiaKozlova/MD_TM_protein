part_analysis <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
paste0(part_analysis)

df_merge<-read.csv(paste0(part_analysis,"din/","df_merge_structure_log.csv"),stringsAsFactors = F)
df_merge<-df_merge%>%select(name.x,receptor,ligand,size_of_group)
df_merge<-unique(df_merge)
if(dir.exists(paste0(part_analysis,"din/complex_structure_surf"))) {system(command = paste0("rm -r ",part_analysis,"din/complex_structure_surf"),ignore.stdout=T,wait = T)}

if(!dir.exists(paste0(part_analysis,"din/complex_structure_surf"))){dir.create(paste0(part_analysis,"din/complex_structure_surf"))}
for (i in 1:nrow(df_merge)) {
  receptor_name<-paste0(part_analysis,"receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_analysis,"din/str_fin/",df_merge$name.x[i])
  pdb_receptor<-read.pdb(receptor_name)
  pdb_ligand<-read.pdb(ligand_name)
  pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
  write.pdb(pdb_complex,paste0(part_analysis,"din/complex_structure_surf/",df_merge$name.x[i]))
}