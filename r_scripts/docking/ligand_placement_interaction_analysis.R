part_start <- commandArgs(trailingOnly=TRUE)
library(dplyr)
library(ggplot2)
library(bio3d)
part_name<-list.files(paste0(part_start,"MD_analysis/"))

#v_seq<-read.csv(paste0(part_start,"start/sequence/",seq_name[1]))

part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")

setwd(part_name)
df_merge<-read.csv(file = paste0(part_name,"din/df_structure.csv"),stringsAsFactors = F)
df_convert<-df_merge%>%select( receptor,system_name,Membrane,Structure)
df_convert<-unique(df_convert)

if (!dir.exists(paste0(part_name,"din/tost/"))) {dir.create(paste0(part_name,"din/tost/"))}
if (!dir.exists(paste0(part_name,"din/tcl/"))) {dir.create(paste0(part_name,"din/tcl/"))}
if (!dir.exists(paste0(part_name,"plot/"))) {dir.create(paste0(part_name,"plot/"))}
if (!dir.exists(paste0(part_name,"din/all_interactions/"))) {dir.create(paste0(part_name,"din/all_interactions/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/fin_data/docking_data/"))) {
  dir.create(paste0(part_start,"MD_analysis/fin_data/docking_data/"))}
if (!dir.exists(paste0(part_name,"din/tost/"))) {
  dir.create(paste0(part_name,"din/tost/"))}
i<-1
#for (i in 1:nrow(df_merge)) {
#  if (!dir.exists(paste0(part_name,"din/tost/",df_merge$receptor[i],"_",df_merge$ligand[i]))) {
#    dir.create(paste0(part_name,"din/tost/",df_merge$receptor[i],"_",df_merge$ligand[i]))}
#  receptor_name<-paste0(part_name,"receptor_start/",df_merge$receptor[i],".pdb")
#  ligand_name<-paste0(part_name,"din/str_fin/",df_merge$name.x[i])
#  pdb_receptor<-read.pdb(receptor_name)
#  pdb_ligand<-read.pdb(ligand_name)
#  pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
  
  #  pdb<-read.pdb(paste0(part_name,"str_fin/",df_merge$name.x[i]))
#  write.pdb(pdb_complex,paste0(part_name,"din/tost/",df_merge$receptor[i],"_",df_merge$ligand[i],"/",df_merge$structure_order[i],"_",df_merge$name.x[i]))
#}

#df_structure<-read.csv(paste0(part_start,"docking/docking_first/din/df_structure.csv"),stringsAsFactors =  F)
df_structure<-df_merge
#df_structure$type[df_structure$type=="surf"]<-"serf"
#df_interaction<-read.csv(paste0(part_name,"din/interaction_",df_structure$type[1],
#                                "/",df_structure$name.x[1],".csv"),
#                         stringsAsFactors = F)
#df_interaction<-df_interaction%>%mutate(name.x=df_structure$name.x[1])
#df_interaction<-left_join(df_interaction,df_structure,by="name.x",
#                          relationship ="many-to-many")
#df_interaction<-df_interaction%>%select()
for (i in 1:nrow(df_merge)) {
  df_interaction<-read.csv(paste0(part_name,"din/interaction_",df_structure$type[i],
                                  "/",df_structure$name.x[i],".csv"),
                           stringsAsFactors = F)
  df_interaction<-df_interaction%>%mutate(name.x=df_structure$name.x[i])
  df_interaction<-left_join(df_interaction,df_structure,by="name.x",
                            relationship ="many-to-many")
  df_interaction<-df_interaction%>%filter(persent_interactions>0)
  write.csv(df_interaction,paste0(part_name,"din/all_interactions/",df_structure$name.x[i],".csv"),row.names = F)
}
i<-1
df_interaction<-read.csv(paste0(part_name,"din/all_interactions/",df_structure$name.x[i],".csv"),stringsAsFactors =  F)
for (i in 2:nrow(df_merge)) {
  df_interaction_add<-read.csv(paste0(part_name,"din/all_interactions/",df_structure$name.x[i],".csv"),stringsAsFactors =  F)
  df_interaction<-rbind(df_interaction,df_interaction_add)
}
write.csv(df_interaction,paste0(part_start,"MD_analysis/fin_data/docking_statistic/path_aminoacid_interactions.csv"),row.names =  F)