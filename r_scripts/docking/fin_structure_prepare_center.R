part_start <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
paste0(part_start)

df_merge<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/din/","df_merge_center.csv"),stringsAsFactors = F)
df_merge<-df_merge%>%mutate(name=paste0(ligand_center,"_",grop_number,"_",models.x))
#df_merge<-semi_join(df_merge,df_all)
df_merge<-df_merge%>%select(name,receptor,ligand,center,size_of_group)
df_merge<-unique(df_merge)
i<-1
if(!dir.exists("complex_structure_center")){dir.create("complex_structure_center")}
for (i in 1:nrow(df_merge)) {
  receptor_name<-paste0(part_start,"MD_analysis/docking/docking_first/receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_start,"MD_analysis/docking/docking_first/din/str_fin/",df_merge$name[i])
  pdb_receptor<-read.pdb(receptor_name)
  pdb_ligand<-read.pdb(ligand_name)
  pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
  write.pdb(pdb_complex,paste0("MD_analysis/docking/docking_first/din/","complex_structure_center/",df_merge$receptor[i],"_",df_merge$ligand[i],"_",
                               df_merge$center[i],"_",df_merge$size_of_group[i],".pdb"))
}
