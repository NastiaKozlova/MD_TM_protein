part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)
v_parta<-list.files("MD")

v_part<-paste0(part_start,"MD/",v_parta)
#part<-v_part[1]
num_model<-1
i<-1
j<-1

if (!dir.exists(paste0("MD_analysis/docking/"))){dir.create(paste0("MD_analysis/docking/"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor"))){dir.create(paste0("MD_analysis/docking/receptor"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor_start"))){dir.create(paste0("MD_analysis/docking/receptor_start"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor"))){dir.create(paste0("MD_analysis/docking/receptor"))}
max_num<-10
main_part<-c(8)
df_main<-data.frame(matrix(ncol=2,nrow = length(main_part)))
colnames(df_main)<-c("number","frames")
df_main$number<-main_part
for (j in 1:length(v_parta)) {
  part<-paste0(part_start,"MD/",v_parta[j])
  if (!dir.exists(paste0(part,"/din/docking/"))){dir.create(paste0(part,"/din/docking/"))}
  for (q in 1:nrow(df_main)) {
    df_main$frames[q]<-length(list.files(paste0(part,"/din/pdb_second/",df_main$number[q],"/")))-1
  }
  df_main<-df_main%>%filter(frames>0)
  df_main<-df_main%>%filter(number==max(number))
  df_ramachadran<-read.csv(paste0("MD_analysis/din/",v_parta[j],"/",df_main[1],"_time_Ramachadran.csv"),stringsAsFactors = F)
  df_ramachadran<-df_ramachadran%>%filter(ramachadran==min(ramachadran))
  for (i in 1:nrow(df_ramachadran)) {
    pdb<-read.pdb(paste0(part,"/din/pdb_second/",df_main$number[1],"/",df_ramachadran$frame_number[i],".pdb"))
    pdb.int<-atom.select(pdb,"water","ions",operator = "OR",inverse=T)
    pdb_fin<-trim.pdb(pdb,pdb.int)
    write.pdb(pdb_fin,paste0("MD_analysis/docking/receptor_start/",v_parta[j],"_",df_ramachadran$frame_number[i],".pdb"))
  }
}