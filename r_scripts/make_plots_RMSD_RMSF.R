part_start <- commandArgs(trailingOnly=TRUE)
library(dplyr)
library(ggplot2)
library(cowplot)
library(bio3d)
library(bio3d)
library(ggplot2)
library(dplyr)
test_10<-seq(from=0,to=1000,by=10)

part<-paste0(part_start,"MD_analysis/")
parta<-paste0(part_start,"MD/")
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)

df_all_systems$fin_name<-as.character(df_all_systems$fin_name)
df_all_systems_productive_run<-df_all_systems%>%filter(Progress=="Productive RUN")
df_all_systems_done<-df_all_systems%>%filter(Progress=="done")
df_all_systems<-rbind(df_all_systems_productive_run,df_all_systems_done)
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
main_part<-c(8)
i<-nrow(df_all_systems)
main<-main_part[2]
if (!dir.exists(paste0(part,"fin_data"))) {dir.create(paste0(part,"fin_data"))}
if (!dir.exists(paste0(part,"fin_data/frame_data"))) {dir.create(paste0(part,"fin_data/frame_data"))}
if (!dir.exists(paste0(part,"fin_data/str_data"))) {dir.create(paste0(part,"fin_data/str_data"))}
if (!dir.exists(paste0(part,"fin_data/str_data_sep"))) {dir.create(paste0(part,"fin_data/str_data_sep"))}
if (!dir.exists(paste0(part,"fin_plots"))) {dir.create(paste0(part,"fin_plots"))}
if (!dir.exists(paste0(part,"fin_plots/test_docking_plots/"))) {dir.create(paste0(part,"fin_plots/test_docking_plots/"))}
if (!dir.exists(paste0(part,"fin_plots/str_plots"))) {dir.create(paste0(part,"fin_plots/str_plots"))}
if (!dir.exists(paste0(part,"fin_plots/str_XYZ"))) {dir.create(paste0(part,"fin_plots/str_XYZ"))}
if (!dir.exists(paste0(part,"fin_plots/frame_plots"))) {dir.create(paste0(part,"fin_plots/frame_plots"))}
if (!dir.exists(paste0(part,"fin_plots/docking_plots"))) {dir.create(paste0(part,"fin_plots/docking_plots"))}
#if (!dir.exists(paste0(part,"fin_plots/lenght_plots"))) {dir.create(paste0(part,"fin_plots/lenght_plots"))}
for (i in 1:nrow(df_all_systems)) {
  for (main in main_part) {
    if (file.exists(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"))){
      df_ramachadran <-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"),stringsAsFactors = F)
      df_ramachadran<-df_ramachadran%>%mutate(frame=number)
      df_Energy<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_",main,".txt"), header=T, na.strings ="", stringsAsFactors= F)
      df_fin<-left_join(df_ramachadran,df_Energy,by=c("frame"="Frame"))
      df_SASA<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/SASA/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
      colnames(df_SASA)<-c("frame","SASA")
      df_fin<-left_join(df_fin,df_SASA,by=c("frame"))
      df_RMSD<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSD/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
      colnames(df_RMSD)<-c("frame","RMSD")
      df_fin<-left_join(df_fin,df_RMSD,by=c("frame"))
      write.csv(df_fin,paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],"_",main,".csv"),row.names = F)
      df_fin<-df_fin%>%mutate(frame=frame/10)
      p_rmsd<-ggplot(data = df_fin)+
        ggtitle(paste0("RMSD"))+
        labs(y = "RMSD (A)", x = "Time (ns)")+
        geom_line(aes(x = frame, y = RMSD))+
        theme_bw()+coord_flip()+
        geom_hline(yintercept = median(df_fin$RMSD,na.rm = T))
      p_sasa<-ggplot(data = df_fin)+
        labs(title=paste("SASA"),
             x = "Time(ns)", y = "SASA (A^2)") +
        geom_line(aes(x = frame,y=SASA))+
        scale_x_continuous(breaks = test_10, labels =  test_10)+
        theme_bw()+coord_flip()+
        geom_hline(yintercept = median(df_fin$SASA,na.rm = T))
      
      p_Elec<-ggplot(data = df_fin)+
        labs(title=paste("Electrostatic energy"),
             x = "Time(ns)", y = "Electrostatic energy (kcal/mol)") +
        geom_line(aes(x = frame,y=Elec))+
        scale_x_continuous(breaks = test_10, labels =  test_10)+
        theme_bw()+coord_flip()+
        geom_hline(yintercept = median(df_fin$Elec))
      p_ramachadran<-ggplot(data = df_fin)+
        labs(title=paste("Ramachadran"), x = "Time(ns)", y = "Ramachadran") +
        geom_line(aes(x = frame,y=ramachadran))+
        scale_x_continuous(breaks = test_10, labels =  test_10)+
        theme_bw()+coord_flip()+
        geom_hline(yintercept = median(df_fin$ramachadran,na.rm = T))
      
      p_Total<-ggplot(data = df_fin)+
        labs(title=paste("Total energy"),
             x = "Time(ns)", y = "Total energy (kcal/mol)") +
        geom_line(aes(x = frame,y=Total))+
        scale_x_continuous(breaks = test_10, labels =  test_10)+
        theme_bw()+coord_flip()+
        geom_hline(yintercept = median(df_fin$Total,na.rm = T))
      p_VdW<-ggplot(data = df_fin)+
        labs(title=paste("VdW"),
             x = "Time(ns)", y = "VdW (kcal/mol)") +
        geom_line(aes(x = frame,y=VdW))+
        scale_x_continuous(breaks = test_10, labels =  test_10)+
        theme_bw()+coord_flip()+
        geom_hline(yintercept = median(df_fin$VdW,na.rm = T))
      p_rmsd_histo<-ggplot(data = df_fin)+
        ggtitle(paste0("RMSD"))+
        labs(x = "RMSD (A)")+
        geom_freqpoly(aes(x = RMSD))+
        theme_bw()+scale_fill_grey()+
        geom_text(x=median(df_RMSD$RMSD), y=20,label=round(median(df_RMSD$RMSD,na.rm = T),digits = 1))+
        geom_vline(xintercept = median(df_RMSD$RMSD,na.rm = T))
      p_sasa_histo<-ggplot(data = df_fin)+
        labs(title="SASA", x = "SASA (A^2)") +
        geom_freqpoly(aes(x = SASA))+
        theme_bw()+
        geom_text(x=median(df_fin$SASA), y=20,label=round(median(df_fin$SASA,na.rm = T),digits = 1))+
        geom_vline(xintercept = median(df_fin$SASA,na.rm = T))
      p_Elec_histo<-ggplot(data = df_fin)+
        labs(title=paste("Electrostatic energy"),
             y = "", x = "Electrostatic energy (kcal/mol)") +
        geom_freqpoly(aes(x = Elec))+
        theme_bw()+
        geom_text(x=median(df_fin$Elec,na.rm = T), y=20,label=round(median(df_fin$Elec,na.rm = T),digits = 1))+
        geom_vline(xintercept = median(df_fin$Elec,na.rm = T))
      p_ramachadran_histo <- ggplot(data = df_fin)+
        labs(title=paste("Ramachadran"),
             y = "", x = "Ramachadran") +
        geom_freqpoly(aes(x = ramachadran),bins=(max(df_fin$ramachadran)-min(df_fin$ramachadran)))+
        theme_bw()+#coord_flip()+
        geom_text(x=median(df_fin$ramachadran,na.rm = T), y=20,label=round(median(df_fin$ramachadran,na.rm = T),digits = 1))+
        geom_vline(xintercept = median(df_fin$ramachadran,na.rm = T))
      p_Total_histo <- ggplot(data = df_fin)+
        labs(title=paste("Total energy"), y = "", x = "Total energy (kcal/mol)") +
        geom_freqpoly(aes(x = Total))+
        theme_bw()+
        geom_vline(xintercept = median(df_fin$Total,na.rm = T))+
        geom_text(x=median(df_fin$Total,na.rm = T), y=20,label=round(median(df_fin$Total,na.rm = T),digits = 1))
      p_VdW_histo<-ggplot(data = df_fin)+
        labs(title=paste("VdW"), x = "VdW (kcal/mol)") +
        geom_freqpoly(aes(x =VdW))+
        theme_bw()+
        geom_text(x=median(df_fin$VdW,na.rm = T), y=20,label=round(median(df_fin$VdW,na.rm = T),digits = 1))+
        geom_vline(xintercept = median(df_fin$VdW,na.rm = T))
      p_all<-plot_grid(p_sasa,       p_rmsd,      p_ramachadran,      p_Total,   
                       p_sasa_histo, p_rmsd_histo,p_ramachadran_histo,p_Total_histo,  
                       nrow=2,
                       rel_heights = c(4.5,1),
                       labels = c("A","B","C","D",
                                  "", "", "", "E"),align="hv")
      ggsave(p_all,   filename = paste0(part,"fin_plots/frame_plots/",df_all_systems$fin_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
    }
    
    if(file.exists(paste0(part_start,"MD_analysis/docking/docking_first/interaction_fin.csv"))){
      df_docking<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/interaction_fin.csv"),stringsAsFactors = F)
      df_docking<-df_docking%>%mutate(receptor=NA)
      for (receptor in 1:nrow(df_docking)) {
        df_docking$receptor[receptor]<-strsplit(df_docking$system[receptor],split = "_")[[1]][1]
      }
      df_docking<-df_docking%>%mutate(receptor=paste0(receptor))
      #    nrow(df_docking)
      df_docking<-df_docking%>%filter(receptor==paste0(df_all_systems$fin_name[i]))
      pdb_name<-df_docking$system[1]
      df_docking<-df_docking%>%select(resno, number_interactions, receptor, center, ligand, grops, grops_number, persent_interactions, aminoacids)
      df_docking<-unique(df_docking)
      pdb<-read.pdb(paste0(part_start,"MD_analysis/docking/receptor_start/",pdb_name,".pdb"))  
      df_seq<-pdb$atom
      df_seq<-df_seq%>%filter(elety=="CA")
      df_seq<-df_seq%>%select(resno,resid, x, y, z)
      df_seq<-left_join(df_seq,df_docking,by="resno")
      
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
      p_conserv<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = z,colour=conservative))+
        theme_bw()#+coord_flip()
      p_ramachadran<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = z,colour=ramachadran))+
        theme_bw()#+coord_flip()
      p_RMSF<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = z,colour=RMSF))+
        theme_bw()
      p_hbonds<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = z,colour=hbonds))+
        theme_bw()
      p_all<-plot_grid(p_conserv,p_ramachadran,       p_hbonds,      p_RMSF,     
                       nrow=4,  labels = c("A","B","C","D"),align="hv")
      ggsave(p_all,   filename = paste0(part,"fin_plots/str_XYZ/",df_all_systems$fin_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
      p_conserv<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = conservative))+
        theme_bw()#+coord_flip()
      p_ramachadran<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = ramachadran))+
        theme_bw()#+coord_flip()
      p_RMSF<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = RMSF))+
        theme_bw()
      p_hbonds<-ggplot(data = df_seq)+
        geom_point(aes(x = resno, y = hbonds))+
        theme_bw()
      p_all<-plot_grid(p_conserv,p_ramachadran,       p_hbonds,      p_RMSF,     
                       nrow=4,  labels = c("A","B","C","D"),align="hv")
      ggsave(p_all,   filename = paste0(part,"fin_plots/str_plots/",df_all_systems$fin_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
      
      write.csv(df_seq,paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),row.names = F)
      df_seq<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors = F)
      df_seq<-df_seq%>%filter(persent_interactions>95)
      df_seq<-df_seq%>%filter(ramachadran==0)
      p_docking<-ggplot(data = df_seq)+
        geom_text(aes(x = resno, y = resid,colour=conservative,label=aminoacids))+
        theme_bw()+ facet_grid(ligand~receptor)
      ggsave(p_docking,   filename = paste0(part,"fin_plots/docking_plots/",df_all_systems$fin_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
      
    }
  }
}
df_seq<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[1],".csv"),stringsAsFactors =  F)
if (nrow(df_all_systems)>1) {
  for (i in 2:nrow(df_all_systems)) {
    df_seq_add<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
    df_seq<-rbind(df_seq,df_seq_add)
    df_seq<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors = F)

  }
}
df_seq<-df_seq%>%filter(persent_interactions==100)
df_seq<-df_seq%>%filter(ramachadran==0)
df_seq<-df_seq%>%filter(conservative>50)
df_seq<-df_seq%>%filter(hbonds>50)
p_docking<-ggplot(data = df_seq)+
  geom_text(aes(x = resno, y = resid,colour=conservative,label=aminoacids))+
  theme_bw()+ facet_grid(ligand~receptor)
ggsave(p_docking,   filename = paste0(part,"fin_plots/docking_fin.png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 


