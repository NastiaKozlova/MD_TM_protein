part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(bio3d)
library(ggplot2)
library(dplyr)

v_main<-c(8)
main_part<-c(8)

v_part<-list.files("MD")
#name<-v_name[1]
#v_namd<-1000
i<-1
main<-main_part[1]
p<-1
v_sep<-seq(from=-30,to=30,by=5)
v_z<-seq(from=-30,to=30,by=1)
df_pore<-data.frame(matrix(ncol = 2,nrow = length(v_z)))
colnames(df_pore)<-c("z","size")
df_pore$z<-v_z
df_pore$size<-0
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/")
  parta<-paste0(part_start,"MD_analysis/fin_plots/docking_plots/")
  #if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
  setwd(paste0(part,"din/pdb_second"))
  if(!dir.exists(paste0(main,"_pore_size/"))){(dir.create(paste0(main,"_pore_size/")))}
  for (main in main_part) {
    if (file.exists(paste0(main,"_hbonds/frame_",0,".csv"))){ 
      frame_number<-length(list.files(path = paste0(main,"_hbonds/")))
      df_topology<-data.frame(matrix(nrow = frame_number, ncol=3))
      colnames(df_topology)<-c("number","frame_number","ramachadran")
      df_topology$number<-0:(frame_number-1)
      df_topology<-df_topology%>%mutate(frame_number=paste0("frame_",number))
      for (i in 1:nrow(df_topology)) {
        df_test<-read.csv(paste0(main,"_hbonds/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
        df_test<-unique(df_test)
        df_test<-df_test%>%group_by(z)%>%mutate(size=n())
        df_pore_size<-ungroup(df_test)
        df_pore_size<-df_pore_size%>%select(z,size)
        df_pore_size<-unique(df_pore_size)
        df_pore_add<-df_pore[!df_pore$z%in%df_pore_size$z,]
        df_pore_size<-rbind(df_pore_size,df_pore_add)
        write.csv(df_pore_size,paste0(main,"_pore_size/",df_topology$frame_number[i],".csv"),row.names = F)
      }
      df_pore_size<-read.csv(paste0(main,"_pore_size/",df_topology$frame_number[1],".csv"),stringsAsFactors = F)
      df_pore_size<-df_pore_size%>%mutate(frame_number=df_topology$number[1])
      for (i in 2:nrow(df_topology)) {
        df_pore_size_add<-read.csv(paste0(main,"_pore_size/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
        df_pore_size_add<-df_pore_size_add%>%mutate(frame_number=df_topology$number[i])
        df_pore_size<-rbind(df_pore_size,df_pore_size_add)
      }
      df_pore_size<-df_pore_size%>%group_by(z)%>%
        mutate(average_size=mean(size))
      df_pore_size<-df_pore_size%>%mutate(colour="bonding")
      df_pore_size$colour[df_pore_size$average_size==0]<-"not bonding"
#      df_pore_size<-df_pore_size%>%filter(average_size)
      p_test<-ggplot(data=df_pore_size)+
        geom_line(aes(x=z,y=average_size))+
        geom_point(aes(x=z,y=average_size,colour = colour))+
        geom_hline(yintercept=0)+
        scale_x_continuous(breaks = v_sep,labels = v_sep)+
        theme_bw()
      ggsave(p_test,filename = paste0(parta,v_part[p],"_average_pore_size.png"), width = 20, height = 12, units = c("cm"), dpi = 1000 ) 
      
      
      
    }
  }
}