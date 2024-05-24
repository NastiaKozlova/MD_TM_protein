part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)
v_parta<-list.files("MD")
v_RMSD<-5
v_part<-paste0(part_start,"MD/",v_parta)
#part<-v_part[1]
num_model<-1
i<-1
j<-3

if (!dir.exists(paste0("MD_analysis/docking/"))){dir.create(paste0("MD_analysis/docking/"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor"))){dir.create(paste0("MD_analysis/docking/receptor"))}
if (!dir.exists(paste0("MD_analysis/docking/receptor_start"))){dir.create(paste0("MD_analysis/docking/receptor_start"))}
#if (!dir.exists(paste0("MD_analysis/docking/receptor_test"))){dir.create(paste0("MD_analysis/docking/receptor_test"))}
#if (!dir.exists(paste0("MD_analysis/docking/df_RMSD"))){dir.create(paste0("MD_analysis/docking/df_RMSD"))}
max_num<-10
main_part<-c(8)
df_main<-data.frame(matrix(ncol=2,nrow = length(main_part)))
colnames(df_main)<-c("number","frames")
df_main$number<-main_part
j<-1
part<-paste0(part_start,"MD/",v_parta[j])
if (!dir.exists(paste0(part,"/din/docking/"))){dir.create(paste0(part,"/din/docking/"))}
pdb_name<-list.files("MD_analysis/cluster_data/claster_pdb/")
for (j in 1:length(pdb_name)) {
  pdb<-read.pdb(paste0("MD_analysis/cluster_data/claster_pdb/",pdb_name[j]))
  write.pdb(pdb,paste0("MD_analysis/docking/receptor_start/",pdb_name[j]))
}
