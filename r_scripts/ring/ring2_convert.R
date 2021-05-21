part_start <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)
test_10<-seq(from=0,to=1000,by=10)

setwd(part_start)
v_part<-list.files("MD")
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/ring2/"))){dir.create(paste0(part_start,"MD_analysis/ring2/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/ring2/script"))){dir.create(paste0(part_start,"MD_analysis/ring2/script"))}
p<-1
i<-1
y<-1
main_part<-c(8)
main<-main_part[1]
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/din/pdb_second/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  setwd(part)
  for (main in main_part) {
    if(file.exists(paste0(parta,main,"_time_Ramachadran.csv"))){
      df_topology<-read.csv(paste0(parta,main,"_time_Ramachadran.csv"),stringsAsFactors = F)
      if (!dir.exists(paste0(main,"_ring/"))){dir.create(paste0(main,"_ring/"))}
      for (i in 1:nrow(df_topology)) {
        if (file.exists(paste0(main,"_ring2/",df_topology$frame_number[i],".txt"))){
          if (!file.exists(paste0(main,"_ring/",df_topology$frame_number[i],".txt"))){
            df_ring<-read_delim(paste0(main,"_ring2/",df_topology$frame_number[i],".txt"), delim = "\t", skip = 11)
            df_ring<-df_ring[1:(which(df_ring$NodeId1%in%"NodeId")-1),]
            df_ring<-df_ring%>%select(NodeId1,Interaction,NodeId2,Distance,Angle,Energy)
            write.csv(df_ring,paste0(main,"_ring/",df_topology$frame_number[i],".txt"),row.names = F)
          }
        }
      }
    }
  }
}
i<-1
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/din/pdb_second/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  setwd(part)
  for (main in main_part) {
    if(file.exists(paste0(parta,main,"_time_Ramachadran.csv"))){
      df_topology<-read.csv(paste0(parta,main,"_time_Ramachadran.csv"),stringsAsFactors = F)
      if (!dir.exists(paste0(main,"_ring_mod/"))){dir.create(paste0(main,"_ring_mod/"))}
      for (i in 1:nrow(df_topology)) {
        if (file.exists(paste0(main,"_ring/",df_topology$frame_number[i],".txt"))){
          if (!file.exists(paste0(main,"_ring_mod/",df_topology$frame_number[i],".txt"))){
            df_ring<-read.csv(paste0(main,"_ring/",df_topology$frame_number[i],".txt"),stringsAsFactors = F)   
            df_ring<-df_ring%>%select(NodeId1,NodeId2,Interaction)
            df_ring<-unique(df_ring)
      
            for (j in 1:nrow(df_ring)) {
              df_ring$NodeId1[j]<-strsplit(df_ring$NodeId1[j],split = ":",fixed = T)[[1]][2]
              df_ring$NodeId2[j]<-strsplit(df_ring$NodeId2[j],split = ":",fixed = T)[[1]][2]
              df_ring$Interaction[j]<-strsplit(df_ring$Interaction[j],split = ":",fixed = T)[[1]][1]
            }
          df_ring<-df_ring%>%mutate(bond=paste(NodeId1,NodeId2,sep = "-"))
          write.csv(df_ring,paste0(main,"_ring_mod/",df_topology$frame_number[i],".txt"),row.names = F)
          }
        }
      }
    }
  }
}