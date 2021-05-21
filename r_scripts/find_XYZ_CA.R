part_start <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)

v_main<-c(8)
setwd(part_start)
v_part<-list.files("MD")
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}
test_10<-seq(from=0,to=1000,by=10)
main_part<-c(8)
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
  setwd(paste0(part,"din/"))
  if (!dir.exists(paste0("pdb_sec/"))){dir.create(paste0("pdb_sec/"))}
  for (main in main_part) {
    number_frame<-length(list.files(paste0("pdb_second/",main)))-1
    if (number_frame>0) { 
      if (!dir.exists(paste0("pdb_sec/",main))){dir.create(paste0("pdb_sec/",main))}
      if (!dir.exists(paste0("pdb_sec/",main,"_SG"))){dir.create(paste0("pdb_sec/",main,"_SG"))}
      for (q in 0:number_frame) {
        pdb<-read.pdb(paste0("pdb_second/",main,"/frame_",q,".pdb"))
        df_pdb<-pdb$atom
        df_pdb<-df_pdb%>%filter(elety=="CA")
        df_pdb<-df_pdb%>%select(resid,resno,x,y,z)
        
        write.csv(df_pdb,paste0("pdb_sec/",main,"/frame_",q,".csv"),row.names = F)
        df_pdb<-pdb$atom
        df_pdb<-df_pdb%>%filter(elety=="SG")
        df_pdb<-df_pdb%>%select(resid,resno,x,y,z)
        write.csv(df_pdb,paste0("pdb_sec/",main,"_SG/frame_",q,".csv"),row.names = F)
      }
    }
  }
}
