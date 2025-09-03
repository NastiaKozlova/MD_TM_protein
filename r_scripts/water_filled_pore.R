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
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/")
  #parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  #if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
  setwd(paste0(part,"din/pdb_second"))
  
  for (main in main_part) {
    if (file.exists(paste0("hbond_",main,"/frame_",0,".pdb"))){ 
      frame_number<-length(list.files(path = paste0(main)))
      df_topology<-data.frame(matrix(nrow = frame_number, ncol=3))
      colnames(df_topology)<-c("number","frame_number","ramachadran")
      df_topology$number<-0:(frame_number-1)
      df_topology<-df_topology%>%mutate(frame_number=paste0("frame_",number))
      if(!dir.exists(paste0(main,"_hbonds/"))){(dir.create(paste0(main,"_hbonds/")))}
      if(!dir.exists(paste0(main,"_water/"))){(dir.create(paste0(main,"_water/")))}
      for (i in 1:nrow(df_topology)) {
        pdb<-read.pdb(paste0("hbond_",main,"/",df_topology$frame_number[i],".pdb"))
        df_pdb<-pdb$atom
        
        df_pdb<-df_pdb%>%mutate(type="membrane")
        df_pdb$type[df_pdb$resid%in%c("MET", "TYR", "LEU", "LYS", "ASN", "THR", "PHE", "TRP", "GLY", "ILE", "ALA", "PRO", 
                                      "HSD", "ASP", "SER", "GLN","ARG","VAL","CYS","GLU")]<-"protein"
        df_pdb$type[df_pdb$resid%in%c("TIP3", "SOD",  "CLA")]<-"water"
        df_pdb<-df_pdb%>%mutate(x=round(x,digits = 0))
        df_pdb<-df_pdb%>%mutate(y=round(y,digits = 0))
        df_pdb<-df_pdb%>%mutate(z=round(z,digits = 0))
        df_pdb<-unique(df_pdb)
        df_water<-df_pdb%>%filter(type=="water")
        write.csv(df_water,paste0(main,"_water/",df_topology$frame_number[i],".csv"),row.names = F)
        df_pdb<-df_pdb%>%group_by(type,x,y,z)%>%mutate(separate_quantity=n())
        df_pdb<-ungroup(df_pdb)
        df_pdb<-df_pdb%>%group_by(x,y,z)%>%mutate(full_quantity=n())
        df_pdb<-ungroup(df_pdb)
        df_pdb<-df_pdb%>%mutate(persent=separate_quantity/full_quantity)%>%
          mutate(max_persent=max(persent))%>%filter(max_persent==persent)
        df_pdb<-df_pdb%>%select(type,x,y,z)%>%mutate(c=NA)
        df_pdb_protein<-df_pdb%>%filter(type=="protein")
        min_x<-round(min(df_pdb_protein$x)-12,digits = 0)
        max_x<-round(max(df_pdb_protein$x)+12,digits = 0)
        
        min_y<-round(min(df_pdb_protein$y)-12,digits = 0)
        max_y<-round(max(df_pdb_protein$y)+12,digits = 0)
        
        min_z<-round(min(df_pdb_protein$z)-12,digits = 0)
        max_z<-round(max(df_pdb_protein$z)+12,digits = 0)
        
        df_pdb<-df_pdb%>%filter(x>min_x)%>%
          filter(y>min_y)%>%
          filter(z>min_z)%>%
          filter(x<max_x)%>%
          filter(y<max_y)%>%
          filter(z<max_z)
        
        #%>%mutate(c=NA)
        df_pdb_water<-df_pdb%>%filter(type=="water")%>%mutate(c=NA)
        df_pdb_merge<-left_join(df_pdb_water,df_pdb_protein,by="c",relationship = "many-to-many")
        df_pdb_merge<-df_pdb_merge%>%mutate(length=sqrt((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2))
        df_pdb_merge<-df_pdb_merge%>%filter(length<12)
        
        df_pdb_merge<-df_pdb_merge%>%select(type.x,x.x,y.x,z.x,c)
        df_pdb_merge<-unique(df_pdb_merge)
        #df_pdb_merge
        colnames(df_pdb_merge)<-colnames(df_pdb)
        #df_pdb<-ungroup(df_pdb)
        max_membrane<-max(df_pdb$z[df_pdb$type=="membrane"])-10
        min_membrane<-min(df_pdb$z[df_pdb$type=="membrane"])+10
        df_test<-df_pdb_merge%>%filter(z<max_membrane)%>%
          filter(z>min_membrane)%>%filter(type=="water")
        write.csv(df_test,paste0(main,"_hbonds/",df_topology$frame_number[i],".csv"),row.names = F)
      }
    }
  }
}
    