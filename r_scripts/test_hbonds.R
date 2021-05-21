part_start <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)

setwd(part_start)
v_part<-list.files("MD")
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}

sort_hbonds<-function(file_name,n_frame){
  df_hbonds<-read_tsv(file_name,skip = 2,col_names = F)
  colnames(df_hbonds)<-c("donor","acceptor","occupancy")
  df_hbonds<-df_hbonds%>%mutate(frame=q)
  df_hbonds<-df_hbonds%>%mutate(amino=NA)
  df_hbonds<-df_hbonds%>%mutate(number=NA)
  df_hbonds<-df_hbonds%>%mutate(type=NA)
  for (i in 1:nrow(df_hbonds)) {
    dd<-strsplit(df_hbonds$donor[i],split = "-",fixed = T)[[1]][1]
    dd<-strsplit(dd,split = "",fixed = T)[[1]]
    number_d<-paste0(dd[4:length(dd)],collapse = "")
    amino_d<-paste0(dd[1:3],collapse = "")
    aa<-strsplit(df_hbonds$acceptor[i],split = "-",fixed = T)[[1]][1]
    aa<-strsplit(aa,split = "",fixed = T)[[1]]
    number_a<-paste0(aa[4:length(aa)],collapse = "")
    amino_a<-paste0(aa[1:3],collapse = "")
    if (amino_a=="TIP") {
      df_hbonds$type[i]<-"donor"
      df_hbonds$number[i]<-number_d
      df_hbonds$amino[i]<-amino_d
    }
    if (amino_d=="TIP") {
      df_hbonds$type[i]<-"acceptor"
      df_hbonds$number[i]<-number_a
      df_hbonds$amino[i]<-amino_a
    }
  }
  df_hbonds<-df_hbonds%>%select(number,amino,number,type,frame)
  df_hbonds<-unique(df_hbonds)
  df_hbonds<-df_hbonds%>%mutate(number=as.numeric(number))
  return(df_hbonds)
}
test_10<-seq(from=0,to=1000,by=10)
main_part<-c(8)
main<-main_part[1]
p<-1
q<-0
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
  setwd(paste0(part,"din/"))

  for (main in main_part) {
    number_frame<-length(list.files(paste0("hbonds/",main,"")))-1
    if (!dir.exists(paste0("hbonds/",main,"_mod"))){dir.create(paste0("hbonds/",main,"_mod"))}
    if (number_frame>0) { 

      for (q in 0:number_frame) {
        file_name<-paste0("hbonds/",main,"/frame_",q,".txt")
        df_bonds<-sort_hbonds(file_name=file_name,n_frame=q)
        df_pdb<-read.csv(paste0("pdb_sec/",main,"/frame_",q,".csv"),stringsAsFactors = F)
        df_bonds<-left_join(df_bonds,df_pdb,c("number"="resno","amino"="resid"))
        df_bonds<-unique(df_bonds)
        write.csv(df_bonds,paste0("hbonds/",main,"_mod/frame_",q,".txt"),row.names = F)
      }
      df_hbonds<-read.csv(paste0("hbonds/",main,"_mod/frame_",0,".txt"),stringsAsFactors = F)
      for (q in 1:number_frame) {
        df_hbonds_add<-read.csv(paste0("hbonds/",main,"_mod/frame_",q,".txt"),stringsAsFactors = F)
        df_hbonds<-rbind(df_hbonds,df_hbonds_add)
      }
      df_hbonds<-df_hbonds%>%filter(abs(z)<18)
      df_hbonds<-df_hbonds%>%select(number, amino,frame)
      df_bonds<-unique(df_hbonds)
      df_bonds<-df_bonds%>%group_by(number)%>%mutate(occupancy=n())
      df_bonds<-df_bonds%>%mutate(persent=occupancy/(number_frame+1)*100)
      df_bonds<-df_bonds%>%select(number, amino, occupancy,	persent)
      df_bonds<-unique(df_bonds)
      write.csv(df_bonds,paste0(part_start,"MD_analysis/din/",v_part[p],"/hbonds_",main,".csv"),row.names = F)                     
    }
  }
}
