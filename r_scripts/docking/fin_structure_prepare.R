part_start <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
paste0(part_start)
df_all<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%select(receptor,ligand)
df_all<-unique(df_all)
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/din/")
setwd(part_name)
i<-1
#for (i in 1:nrow(df_all)) {
#  if(!file.exists(paste0("interaction_serf/",df_all$receptor[i],"_",df_all$ligand[i],".csv"))){
#    df_all$receptor[i]<-NA
#  }
#}
#df_all<-df_all%>%filter(!is.na(receptor))

v_complex<-list.files(paste0("interaction_serf/"))
a<-c()
i<-1
for (i in 1:length(v_complex)) {
  b<-strsplit(v_complex[i],split = ".",fixed = T)[[1]][1]
  b<-paste0(b,".pdb")
  a<-c(a,b)
}
df_merge<-read.csv(paste0(part_name,"df_merge_structure_log.csv"),stringsAsFactors = F)
df_merge<-df_merge[df_merge$name.x%in%a,]
#df_merge<-semi_join(df_merge,df_all)
df_merge<-df_merge%>%select(name.x,receptor,ligand, size_of_group)
df_merge<-unique(df_merge)
i<-1
if(!dir.exists("complex_structure")){dir.create("complex_structure")}
for (i in 1:nrow(df_merge)) {
  receptor_name<-paste0(part_start,"MD_analysis/docking/docking_first/receptor_start/",df_merge$receptor[i],".pdb")
  if(file.exists(paste0(part_name,"structure_merged/",df_merge$name.x[i]))){
    ligand_name<-paste0(part_name,"structure_merged/",df_merge$name.x[i])
    pdb_receptor<-read.pdb(receptor_name)
    pdb_ligand<-read.pdb(ligand_name)
    pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
    write.pdb(pdb_complex,paste0("complex_structure/",df_merge$name.x[i]))
  }
}
