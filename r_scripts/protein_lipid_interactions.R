#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
#install.packages("stringr")
library(stringr)
library(bio3d)
library(dplyr)
setwd(part_start)
v_parta<-list.files('MD')
v_part<-paste0(part_start,'MD/',v_parta)
#namd_exe<-paste0(part_start,"programs/NAMD_2.14_Linux-x86_64-multicore/namd2 ")

#if (!dir.exists(paste0(part_start,'MD_analysis/tcl/'))){dir.create(paste0(part_start,'MD_analysis/tcl/'))}
#if (!dir.exists(paste0(part_start,'MD_analysis/hbonds/'))){dir.create(paste0(part_start,'MD_analysis/hbonds/'))}
if (!dir.exists(paste0(part_start,'MD_analysis/din/'))){dir.create(paste0(part_start,'MD_analysis/din/'))}
i<-1
q<-1
j<-4
for (j in 1:length(v_parta)) {
  part<-paste0(part_start,'MD/',v_parta[j],"/")
  #    if(file.exists(paste0(part,"namd/step",8,".dcd"))){
  #        parta<-paste0(v_parta)
  
  #        if (!dir.exists(paste0(part,'/din'))){dir.create(paste0(part,'/din'))}
  #        if (!dir.exists(paste0(part,'/din/pdb_second'))){dir.create(paste0(part,'/din/pdb_second'))}
  #        if (!dir.exists(paste0(part,'/din/pdb_second/',8))){dir.create(paste0(part,'/din/pdb_second/',8))}
  #        if (!dir.exists(paste0(part,'/din/pdb_second/hbond_',8))){dir.create(paste0(part,'/din/pdb_second/hbond_',8))}
  
  #        if (!dir.exists(paste0(part,'/din/RMSD'))){dir.create(paste0(part,'/din/RMSD'))}
  #        if (!dir.exists(paste0(part,'/din/RMSF'))){dir.create(paste0(part,'/din/RMSF'))}
  #        if (!dir.exists(paste0(part,'/din/SASA'))){dir.create(paste0(part,'/din/SASA'))}
  #        if (!dir.exists(paste0(part,'/din/Energy'))){dir.create(paste0(part,'/din/Energy'))}
  #        if (!dir.exists(paste0(part,'/din/Energy_protein_lipid'))){dir.create(paste0(part,'/din/Energy_protein_lipid'))}
  #        if (!dir.exists(paste0(part,'/din/hbonds'))){dir.create(paste0(part,'/din/hbonds'))}
  if (!dir.exists(paste0(part_start,'MD_analysis/din/',v_parta[j],"/"))){
    dir.create(paste0(part_start,'MD_analysis/din/',v_parta[j],"/"))}
  pdb<-read.pdb(paste0(part,"din/pdb_second/8/frame_0.pdb"))
  #        pdb.int<-atom.select(pdb,"protein", inverse=TRUE)
  #        pdb <- trim.pdb(pdb, pdb.int)
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-df_pdb%>%select(resno,segid,x,y,z)
  df_pdb<-unique(df_pdb)
  q<-1
  df_energy<-read.table(paste0(part,'din/Energy_protein_lipid/',df_pdb$resno[q],
                               '_segname_',df_pdb$segid[q],'.txt'), header=T,
                        na.strings ="", stringsAsFactors= F)
  df_energy<-df_energy%>%mutate(resno=df_pdb$resno[q])
  df_energy<-df_energy%>%mutate(segid=df_pdb$segid[q])
  if(nrow(df_pdb)>1){
    for (q in 2:nrow(df_pdb)) {
      df_energy_add<-read.table(paste0(part,'din/Energy_protein_lipid/',df_pdb$resno[q],
                                       '_segname_',df_pdb$segid[q],'.txt'), header=T,
                                na.strings ="", stringsAsFactors= F)
      df_energy_add<-df_energy_add%>%mutate(resno=df_pdb$resno[q])
      df_energy_add<-df_energy_add%>%mutate(segid=df_pdb$segid[q])
      df_energy<-rbind(df_energy,df_energy_add)
      
    }
  }
  
  
  df_energy<-left_join(df_energy,df_pdb,by = c("resno", "segid"))
  df_energy<-unique(df_energy)
  write.csv(df_energy,paste0(part_start,"MD_analysis/din/",v_parta[j],"/Energy_protein_lipid.csv"),row.names = F)
  
  df_energy<-df_energy%>%group_by(resno,segid)%>%mutate(number=n())
  df_energy<-df_energy%>%select(resno,segid,number)
  print(unique(df_energy$number))
  
}
