part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(dplyr)
library(Peptides)
library(ggplot2)
library(cowplot)
library(bio3d)

test_10<-seq(from=0,to=1000,by=50)
v_pallete<-c("sheet"="#BBBBBB","helix"="#333333")
df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors =  F)
part<-paste0(part_start,"MD_analysis/")
parta<-paste0(part_start,"MD/")
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)

df_all_systems$system_name<-as.character(df_all_systems$system_name)
#df_all_systems<-df_all_systems%>%mutate(Progress=="notdone")
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
main_part<-c(8)
i<-nrow(df_all_systems)
main<-main_part[2]


if (dir.exists(paste0(part,"fin_data"))) { system(command = paste0("rm -r ",part,"fin_data"))}
if (dir.exists(paste0(part,"fin_plots"))) { system(command = paste0("rm -r ",part,"fin_plots"))}
if (!dir.exists(paste0(part,"fin_data"))) {dir.create(paste0(part,"fin_data"))}
if (!dir.exists(paste0(part,"fin_data/frame_data"))) {dir.create(paste0(part,"fin_data/frame_data"))}
if (!dir.exists(paste0(part,"fin_data/frame_statisitc"))) {dir.create(paste0(part,"fin_data/frame_statisitc"))}
if (!dir.exists(paste0(part,"fin_data/str_data"))) {dir.create(paste0(part,"fin_data/str_data"))}
if (!dir.exists(paste0(part,"fin_data/aminoacids_interactions/"))){dir.create(paste0(part,"fin_data/aminoacids_interactions/"))}

if (!dir.exists(paste0(part,"fin_plots"))) {dir.create(paste0(part,"fin_plots"))}
if (!dir.exists(paste0(part,"fin_plots/str_plots"))) {dir.create(paste0(part,"fin_plots/str_plots"))}
if (!dir.exists(paste0(part,"fin_plots/frame_plots"))) {dir.create(paste0(part,"fin_plots/frame_plots"))}
if (!dir.exists(paste0(part,"fin_plots/frame_statisitc"))) {dir.create(paste0(part,"fin_plots/frame_statisitc"))}
if (!dir.exists(paste0(part,"fin_plots/docking_plots"))) {dir.create(paste0(part,"fin_plots/docking_plots"))}

i<-1
main<-main_part[1]
for (i in 1:nrow(df_all_systems)) {
  for (main in main_part) {
    if (file.exists(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_",main,".txt"))){
      if (file.exists(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"))){
        if (file.exists(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSF/",main,".txt"))){
          df_all_systems$Progress[i]<-"done"
#          df_second_structure<-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/Second_structure_",main,".csv"),stringsAsFactors = F)
#          df_second_structure<-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/Second_structure_",main,".csv"),stringsAsFactors = F)
          
          df_RMSF <- read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSF/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
          colnames(df_RMSF)<-"RMSF"
          df_RMSF <- df_RMSF%>% mutate(Resid=c(1:nrow(df_RMSF)))
          
          df_ramachadran <-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"),stringsAsFactors = F)
          df_ramachadran<-df_ramachadran%>%mutate(frame=number)
          df_Energy_protein<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_",main,".txt"), header=T, na.strings ="", stringsAsFactors= F)
          colnames(df_Energy_protein)<-c(colnames(df_Energy_protein)[1:2],paste0("protein_",colnames(df_Energy_protein)[3:ncol(df_Energy_protein)]))
          df_Energy_protein_water<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_water_energy_",main,".txt"), header=T, na.strings ="", stringsAsFactors= F)
          colnames(df_Energy_protein_water)<-c(colnames(df_Energy_protein_water)[1:2],paste0("protein_water_",colnames(df_Energy_protein_water)[3:ncol(df_Energy_protein_water)]))
          df_Energy_protein_lipids<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_lipid_energy_",main,".txt"), header=T, na.strings ="", stringsAsFactors= F)
          colnames(df_Energy_protein_lipids)<-c(colnames(df_Energy_protein_lipids)[1:2],paste0("protein_lipid_",colnames(df_Energy_protein_lipids)[3:ncol(df_Energy_protein_lipids)]))
          df_frame_data<-full_join(df_ramachadran,df_Energy_protein,by=c("frame"="Frame"))
          df_frame_data$time<-NULL
          df_frame_data$Time<-NULL
          df_frame_data<-left_join(df_frame_data,df_Energy_protein_water,by=c("frame"="Frame"))
          df_frame_data$Time<-NULL
          df_frame_data$time<-NULL
          df_frame_data<-left_join(df_frame_data,df_Energy_protein_lipids,by=c("frame"="Frame"))
          df_frame_data$Time<-NULL
          df_frame_data$time<-NULL
          df_SASA<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/SASA/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
          colnames(df_SASA)<-c("frame","SASA_protein")
          df_frame_data<-left_join(df_frame_data,df_SASA,by=c("frame"))
          df_frame_data$Time<-NULL
          df_frame_data$time<-NULL
          df_RMSD_protein<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSD/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
          colnames(df_RMSD_protein)<-c("frame","RMSD_protein")
          df_frame_data<-left_join(df_frame_data,df_RMSD_protein,by=c("frame"))
          df_frame_data$Time<-NULL
          df_frame_data$time<-NULL 
          df_frame_data<-df_frame_data%>%mutate(frame=frame/10)
          write.csv(df_frame_data,paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],".csv"),row.names = F)

        }
      }
    }
  }
}
i<-1
for (i in 1:nrow(df_all_systems)) {
  for (main in main_part) {
    if (file.exists(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"))){
      df_ramachadran_time <-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"),stringsAsFactors = F)
      if(file.exists(paste0(part_start,"MD_analysis/docking/receptor_start/",df_all_systems$fin_name[i],".pdb"))){
        pdb<-read.pdb(paste0(part_start,"MD_analysis/docking/receptor_start/",df_all_systems$fin_name[i],".pdb")) 
        df_seq<-pdb$atom
        df_seq<-df_seq%>%filter(elety=="CA")
        df_seq<-df_seq%>%select(resno,resid, x, y, z,type)
        for (q in 1:nrow(df_topology)) {
          df_seq$type[(df_seq$resno>=df_topology$seq_beg[q])&(df_seq$resno<=df_topology$seq_end[q])]<-df_topology$type[q]
          df_seq$topology[(df_seq$resno>=df_topology$seq_beg[q])&(df_seq$resno<=df_topology$seq_end[q])]<-df_topology$topology[q]
        }
        
        df_seq<-df_seq%>%mutate(ramachadran=0)
        df_ramachadran <-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_ramachad.csv"),stringsAsFactors = F)
        for (p in 1:nrow(df_ramachadran)) {
          df_seq$ramachadran[df_seq$resno==df_ramachadran$number[p]]<-df_ramachadran$x_num[p]
        }
        df_hbonds <-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/hbonds_",main,".csv"),stringsAsFactors = F)
        df_seq<-df_seq%>%mutate(hbonds=0)
        for (p in 1:nrow(df_hbonds)) {
          df_seq$hbonds[df_seq$resno==df_hbonds$number[p]]<-df_hbonds$persent[p]
        }
        df_RMSF<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSF/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
        colnames(df_RMSF)<-"RMSF"
        df_RMSF<-df_RMSF%>%mutate(resno=1:nrow(df_RMSF))
        df_seq<-left_join(df_seq,df_RMSF,by=c("resno"))
        v_conserv<-list.files(paste0(part_start,"MD_analysis/conservative/"))
        df_conserv<-read.csv(paste0(part_start,"MD_analysis/conservative/",v_conserv),stringsAsFactors =  F)
        df_seq<-left_join(df_seq,df_conserv,by=c("resno"="number"))
        
        df_ring<-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/8_ring_interaction.csv"),stringsAsFactors = F)
        df_ring<-df_ring%>%mutate(persent_intractions=all/nrow(df_ramachadran_time))
        df_ring<-df_ring%>%mutate(persent_HBOND=HBOND/nrow(df_ramachadran_time))
        df_ring<-df_ring%>%mutate(persent_VDW=VDW/nrow(df_ramachadran_time))
        df_ring<-df_ring%>%mutate(persent_PIPISTACK=PIPISTACK/nrow(df_ramachadran_time))
        df_ring<-df_ring%>%mutate(persent_IONIC=IONIC/nrow(df_ramachadran_time))
        df_ring<-df_ring%>%mutate(persent_IAC=IAC/nrow(df_ramachadran_time))
        df_ring<-df_ring%>%mutate(persent_PICATION=PICATION/nrow(df_ramachadran_time))
        df_ring1<-df_ring%>%select(number.x,number.y,persent_intractions,
                                   persent_HBOND,persent_VDW,persent_PIPISTACK,
                                   persent_IONIC,persent_IAC,persent_PICATION)
        df_ring2<-df_ring%>%select(number.y,number.x,persent_intractions,
                                   persent_HBOND,persent_VDW,persent_PIPISTACK,
                                   persent_IONIC,persent_IAC,persent_PICATION)
        colnames(df_ring2)<-colnames(df_ring1)
        df_ring<-rbind(df_ring1,df_ring2)
        write.csv(df_ring,paste0(part,"fin_data/aminoacids_interactions/",df_all_systems$fin_name[i],".csv"),row.names = F)
        
        df_ring<-df_ring%>%filter(persent_intractions>0)
        df_ring<-df_ring%>%filter(persent_intractions>quantile(df_ring$persent_intractions,probs = 0.75))
        df_ring<-df_ring%>%filter(abs(number.y-number.x)>5)
        df_ring<-df_ring%>%group_by(number.x)%>%mutate(importance_interprotein_intractions=n())
        
        
        df_ring<-df_ring%>%select(number.x,importance_interprotein_intractions)
        df_ring<-unique(df_ring)
        
        df_seq<-left_join(df_seq,df_ring,by=c("resno"="number.x"))
        df_seq$importance_interprotein_intractions[is.na(df_seq$importance_interprotein_intractions)]<-0
        df_seq<-df_seq%>%mutate(charge=charge(amino,pH = 7.4,pKscale = "Stryer"))
        df_seq<-df_seq%>%mutate(hydrophobicity=hydrophobicity(amino,scale = "Engelman"))
        write.csv(df_seq,paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),row.names = F)
      }
    }
  }
}
