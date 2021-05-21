#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
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
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  setwd(part)
  for (y in 1:length(main_part)) {
    main<-main_part[y]
    frame_number<-length(list.files(path = paste0("din/pdb_second/",main,"_rama")))
    if (frame_number>0) {
      df_topology<-data.frame(matrix(nrow = frame_number, ncol=4))
      colnames(df_topology)<-c("number","frame_number","main","system")
      df_topology$system<-v_part[p]
      df_topology$main<-main
      df_topology$number<-0:(frame_number-1)
      df_topology<-df_topology%>%mutate(frame_number=paste0("frame_",number))
      df_topology<-df_topology%>%mutate(script=NA)
      if (!dir.exists(paste0(part,"din/pdb_second/",main,"_ring2/"))){dir.create(paste0(part,"din/pdb_second/",main,"_ring2/"))}
    
      for (i in 1:nrow(df_topology)) {
        df_topology$script[i]<-paste0(part_start,"programs/dist/bin/Ring -i ",part,"din/pdb_second/",main,"/",df_topology$frame_number[i],".pdb > ",part,"din/pdb_second/",main,"_ring2/",df_topology$frame_number[i],".txt\n",
                                      "rm ", part,"din/pdb_second/",main,"/",df_topology$frame_number[i],".pdb_fasta_P\n",
                                      "rm ", part,"din/pdb_second/",main,"/",df_topology$frame_number[i],".pdb_modified\n")
      }
      write.csv(df_topology,paste0(part_start,"MD_analysis/ring2/script/",v_part[p],"_",main,".txt"), row.names = F)
    }
  }
}
v_topology<-list.files(paste0(part_start,"MD_analysis/ring2/script/"))

df_script<-data.frame(matrix(ncol=1,nrow = length(v_topology)))
colnames(df_script)<-c("script")
for (i in 1:length(v_topology)) {
  df_topology<-read.csv(paste0(part_start,"MD_analysis/ring2/script/",v_topology[i]),stringsAsFactors = F)
  df_script$script[i]<-paste0(df_topology$script,collapse = "\n")
}
df_script_add<-data.frame(matrix(ncol=ncol(df_script),nrow = 1))
colnames(df_script_add)<-colnames(df_script)

df_script_add$script[1]<-paste0("cd ", part_start,"programs/dist/")
df_script<-rbind(df_script_add,df_script)
df_script[is.na(df_script)]<-""
write.table(df_script,paste0(part_start,"MD_analysis/ring2/script.txt"), quote=F,row.names = F,col.names = F,sep = "\n")
system(paste0("chmod +x ",part_start,"MD_analysis/ring2/script.txt"),ignore.stdout=T,wait = T)
system(paste0(part_start,"MD_analysis/ring2/script.txt"),ignore.stdout=T,wait = T)
