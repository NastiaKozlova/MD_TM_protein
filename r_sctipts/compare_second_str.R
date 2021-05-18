#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(cowplot)
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}
number_y<-c(0:10000)
v_pallete<-c("sheet"="#BBBBBB","helix"="#333333")
#v_pallete<-c("sheet"="#BBBBBB","helix"="#FF1FFF")

Extract_Secondary_Structure_From_Pdb<-function(pdb, number_of_pdb){
  m<-dssp(pdb)
  df_second_helix<-data.frame(matrix(nrow = length(m$helix$start),ncol = 3))
  df_second_sheet<-data.frame(matrix(nrow = length(m$sheet$start),ncol = 3))
  if(length(m$sheet$start)>0){
    df_second_sheet$X3<-"sheet"
    for (i in 1:length(m$sheet$start)) {
      df_second_sheet$X1<-m$sheet$start
      df_second_sheet$X2<-m$sheet$end
    }
    df_second_sheet<-df_second_sheet%>%mutate(level_min=number_of_pdb-0.5)
    df_second_sheet<-df_second_sheet%>%mutate(level_max=number_of_pdb+0.5)
  }
  if(length(m$helix$start)>0){
    df_second_helix$X3<-"helix"
    for (i in 1:length(m$helix$start)) {
      df_second_helix$X1<-m$helix$start
      df_second_helix$X2<-m$helix$end
    }
    df_second_helix<-df_second_helix%>%mutate(level_min=number_of_pdb-0.5)
    df_second_helix<-df_second_helix%>%mutate(level_max=number_of_pdb+0.5)
  }
  df_second_all<-rbind(df_second_helix,df_second_sheet)
  colnames(df_second_all)<-c("start","finish","type","level_min","level_max")
  return(df_second_all)
}
setwd(part_start)
v_part<-list.files("MD")
main_part<-c("8")
part<-v_part[1]
num_model<-1
for (p in 1:length(v_part)) {
  for (main in main_part) {
  if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
  part<-paste0(part_start,"MD/",v_part[p],"/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  setwd(paste0(part,"/din/pdb_second"))
  if(dir.exists(main)){
    frame_number<-length(list.files(path = main))-1
    if (file.exists(paste0(main,"/frame_",0,".pdb"))){
      pdb<-read.pdb(paste0(main,"/frame_",0,".pdb"))
      protein.inds <- atom.select(pdb, "protein",resno=2:689)
      backpdb <- trim.pdb(pdb, protein.inds)
      df_pdb_all<-Extract_Secondary_Structure_From_Pdb(pdb = backpdb,number_of_pdb = 0)
      for (i in 1:frame_number) {
        pdb<-read.pdb(paste0(main,"/frame_",i,".pdb"))
        protein.inds <- atom.select(pdb, "protein")
        backpdb <- trim.pdb(pdb, protein.inds)
        df_pdb<-Extract_Secondary_Structure_From_Pdb(pdb = backpdb,number_of_pdb = i)
        df_pdb_all<-rbind(df_pdb_all,df_pdb)
        rm(df_pdb)
      }
      df_pdb_all$level_min<-df_pdb_all$level_min/10
      df_pdb_all$level_max<-df_pdb_all$level_max/10
      write.csv(df_pdb_all,file = paste0(parta,"Second_structure_",main,".csv"),row.names = F)
      rm(df_pdb_all)
      df_pdb_all<-read.csv( paste0(parta,"Second_structure_",main,".csv"),stringsAsFactors =  F)
      test_10<-seq(from=0,to=690,by=10)
      p_second<-ggplot(data = df_pdb_all)+
        ggtitle(paste0("Second structure ",main))+
        labs(x = "Number of aminoasids", y = "Time (ns)")+
        geom_rect(aes(xmin = start, ymin = level_min, xmax= finish, ymax = level_max, colour = type,fill=type))+
        scale_color_manual(values = v_pallete)+ scale_fill_manual(values = v_pallete)+
        scale_x_continuous(breaks = test_10, labels =  test_10)+
        scale_y_continuous(breaks = number_y, minor_breaks = NULL)+
        theme_bw()+theme(legend.position = "top")
      ggsave(p_second,filename = paste0(parta,"Second_strucure_evolution_",main,".png"), width = 40, height = 50, units = c("cm"), dpi = 200 ) 
      }
    }
  }
}
