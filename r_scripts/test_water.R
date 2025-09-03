part_start <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)

setwd(part_start)
v_part<-list.files("MD")
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}

sort_water<-function(file_name,n_frame){
  df_water<-read_tsv(file_name,skip = 2,col_names = F)
  colnames(df_water)<-c("donor","acceptor","occupancy")
  df_water<-df_water%>%mutate(frame=n_frame)
  df_water<-df_water%>%mutate(amino=NA)
  df_water<-df_water%>%mutate(number=NA)
  df_water<-df_water%>%mutate(type=NA)
  for (i in 1:nrow(df_water)) {
    dd<-strsplit(df_water$donor[i],split = "-",fixed = T)[[1]][1]
    dd<-strsplit(dd,split = "",fixed = T)[[1]]
    number_d<-paste0(dd[4:length(dd)],collapse = "")
    amino_d<-paste0(dd[1:3],collapse = "")
    aa<-strsplit(df_water$acceptor[i],split = "-",fixed = T)[[1]][1]
    aa<-strsplit(aa,split = "",fixed = T)[[1]]
    number_a<-paste0(aa[4:length(aa)],collapse = "")
    amino_a<-paste0(aa[1:3],collapse = "")
    if (amino_a!="TIP") {
      df_water$type[i]<-"donor"
      df_water$number[i]<-number_a
      df_water$amino[i]<-amino_a
    }
    if (amino_d!="TIP") {
      df_water$type[i]<-"acceptor"
      df_water$number[i]<-number_d
      df_water$amino[i]<-amino_d
    }
  }
  df_water<-df_water%>%select(number,amino,number,type,frame)
  df_water<-unique(df_water)
  df_water<-df_water%>%mutate(number=as.numeric(number))
  return(df_water)
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
    if (!dir.exists(paste0("water/"))){dir.create(paste0("water/"))}
    if (!dir.exists(paste0("water/",main,"_mod"))){dir.create(paste0("water/",main,"_mod"))}
    if (number_frame>0) { 
      
      for (q in 0:number_frame) {
        file_name<-paste0("hbonds/",main,"/frame_",q,".txt")
        df_bonds<-sort_water(file_name=file_name,n_frame=q)
        df_pdb<-read.csv(paste0("pdb_sec/",main,"/frame_",q,".csv"),stringsAsFactors = F)
        df_bonds<-left_join(df_bonds,df_pdb,c("number"="resno","amino"="resid"))
        df_bonds<-unique(df_bonds)
        write.csv(df_bonds,paste0("water/",main,"_mod/frame_",q,".txt"),row.names = F)
      }
 #     df_water<-read.csv(paste0("water/",main,"_mod/frame_",0,".txt"),stringsAsFactors = F)
#      for (q in 1:number_frame) {
#        df_water_add<-read.csv(paste0("water/",main,"_mod/frame_",q,".txt"),stringsAsFactors = F)
#        df_water<-rbind(df_water,df_water_add)
#      }
#      #      df_water<-df_water%>%filter(abs(z)<18)
#      df_water<-df_water%>%select(number, amino,frame)
#      df_bonds<-unique(df_water)
#      df_bonds<-df_bonds%>%group_by(number)%>%mutate(occupancy=n())
#      df_bonds<-df_bonds%>%mutate(persent=occupancy/(number_frame+1)*100)
#      df_bonds<-df_bonds%>%select(number, amino, occupancy,	persent)
#      df_bonds<-unique(df_bonds)
#      write.csv(df_bonds,paste0(part_start,"MD_analysis/din/",v_part[p],"/water_",main,".csv"),row.names = F)                     
    }
  }
}
