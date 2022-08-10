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

setwd(part_start)
v_part<-list.files("MD")
main_part<-c("8")
part<-v_part[1]
main<-main_part[1]
num_model<-1
for (p in 1:length(v_part)) {
  for (main in main_part) {
    if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
    part<-paste0(part_start,"MD/",v_part[p],"/")
    parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
    setwd(paste0(part,"/din/pdb_second"))
    frame_number<-length(list.files(path = main))
      if (file.exists(paste0(main,"/frame_",0,".pdb"))){
        pdb<-read.pdb(paste0(main,"/frame_",0,".pdb"))
        df_pdb<-pdb$atom
        df_pdb<-df_pdb%>%filter(elety=="CA")
        df_pdb<-df_pdb%>%select(type,resid,resno)
        df_pdb<-df_pdb%>%mutate(second_structure=0)
        df_pdb_loop<-df_pdb%>%mutate(type="loop")
        df_pdb_alpha<-df_pdb%>%mutate(type="helix")
        df_pdb_beta<-df_pdb%>%mutate(type="sheet")

        df_pdb_all<-read.csv(file = paste0(parta,"Second_structure_",main,".csv"),stringsAsFactors = F)


#        write.csv(df_pdb_all,file = paste0(parta,"Second_structure_",main,".csv"),row.names = F)
#        rm(df_pdb_all)
        df_pdb_all<-read.csv( paste0(parta,"Second_structure_",main,".csv"),stringsAsFactors =  F)
        df_pdb_all_alpha<-df_pdb_all%>%filter(type=="helix")
        df_pdb_all_beta<-df_pdb_all%>%filter(type=="sheet")
        
        for (i in 1:nrow(df_pdb__all_alpha)) {
          v_str<-c(df_pdb_all_alpha$start[i]:df_pdb__all_alpha$finish[i])
          df_pdb_alpha$second_structure[df_pdb_alpha$resno%in%v_str]<-df_pdb_alpha$second_structure[df_pdb_alpha$resno%in%v_str]+1
        }
        
        for (i in 1:nrow(df_pdb_beta)) {
          v_str<-c(df_pdb_all_beta$start[i]:df_pdb_all_beta$finish[i])
          df_pdb_beta$second_structure[df_pdb_beta$resno%in%v_str]<-df_pdb_beta$second_structure[df_pdb_beta$resno%in%v_str]+1
        }
#        df_pdb<-rbind(df_pdb_alpha,df_pdb_beta)
        for (t in 1:length(df_pdb$resno)) {
          df_pdb_loop$second_structure[df_pdb_loop$resno==df_pdb$resno[t]]<-frame_number-
            df_pdb_alpha$second_structure[df_pdb_alpha$resno==df_pdb$resno[t]]-
            df_pdb_beta$second_structure[df_pdb_beta$resno==df_pdb$resno[t]]
        }
        df_pdb<-rbind(df_pdb_alpha,df_pdb_beta,df_pdb_loop)
        df_pdb<-df_pdb%>%mutate(persernt_second_structure=second_structure/frame_number*100)
        df_pdb<-df_pdb%>%filter(persernt_second_structure>0)
        df_pdb<-df_pdb%>%group_by(resno)%>%mutate(max_persernt_second_structure=max(persernt_second_structure))
        df_pdb<-ungroup(df_pdb)
        df_pdb<-df_pdb%>%filter(max_persernt_second_structure==persernt_second_structure)
#        df_pdb<-df_pdb%>%filter(persernt_second_structure>50)
        test_10<-seq(from=0,to=690,by=10)
        p_second<-ggplot(data = df_pdb)+
          ggtitle(paste0("Second structure ",main))+
          labs(x = "Number of aminoasids", y = "Time (ns)")+
          geom_point(aes(x = resno, y=persernt_second_structure, colour = type))+
#          scale_color_manual(values = v_pallete)+ scale_fill_manual(values = v_pallete)+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
#          scale_y_continuous(breaks = number_y, minor_breaks = NULL)+
          theme_bw()+theme(legend.position = "top")
        ggsave(p_second,filename = paste0(parta,"Second_strucure_evolution_",main,".png"), width = 40, height = 50, units = c("cm"), dpi = 200 ) 
      }
    }
  }
}