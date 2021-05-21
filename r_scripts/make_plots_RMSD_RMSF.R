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
df_all_systems<-read.csv(paste0(part,"all_systems.csv"),stringsAsFactors = F)
df_all_systems_productive_run<-df_all_systems%>%filter(Progress=="Productive RUN")
df_all_systems_done<-df_all_systems%>%filter(Progress=="done")
df_all_systems<-rbind(df_all_systems_productive_run,df_all_systems_done)
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
main_part<-c(6.6,7)
i<-nrow(df_all_systems)
main<-main_part[2]
if (!dir.exists(paste0(part,"fin_data"))) {dir.create(paste0(part,"fin_data"))}
if (!dir.exists(paste0(part,"fin_data/frame_data"))) {dir.create(paste0(part,"fin_data/frame_data"))}
if (!dir.exists(paste0(part,"fin_data/str_data"))) {dir.create(paste0(part,"fin_data/str_data"))}
if (!dir.exists(paste0(part,"fin_data/str_data_sep"))) {dir.create(paste0(part,"fin_data/str_data_sep"))}
if (!dir.exists(paste0(part,"fin_plots"))) {dir.create(paste0(part,"fin_plots"))}
if (!dir.exists(paste0(part,"fin_plots/test_docking_plots/"))) {dir.create(paste0(part,"fin_plots/test_docking_plots/"))}
if (!dir.exists(paste0(part,"fin_plots/str_plots"))) {dir.create(paste0(part,"fin_plots/str_plots"))}
if (!dir.exists(paste0(part,"fin_plots/frame_plots"))) {dir.create(paste0(part,"fin_plots/frame_plots"))}
if (!dir.exists(paste0(part,"fin_plots/docking_plots"))) {dir.create(paste0(part,"fin_plots/docking_plots"))}
#if (!dir.exists(paste0(part,"fin_plots/lenght_plots"))) {dir.create(paste0(part,"fin_plots/lenght_plots"))}
for (i in 1:nrow(df_all_systems)) {
  for (main in main_part) {
    if (file.exists(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"))){
      df_ramachadran <-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"),stringsAsFactors = F)
      df_ramachadran<-df_ramachadran%>%mutate(frame=number)
      df_Energy<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_",main,"_1.txt"), header=T, na.strings ="", stringsAsFactors= F)
      df_fin<-left_join(df_ramachadran,df_Energy,by=c("frame"="Frame"))
      df_SASA<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/SASA/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
      colnames(df_SASA)<-c("frame","SASA")
      df_fin<-left_join(df_fin,df_SASA,by=c("frame"))
      df_RMSD<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSD/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
      colnames(df_RMSD)<-c("frame","RMSD")
      df_fin<-left_join(df_fin,df_RMSD,by=c("frame"))
      write.csv(df_fin,paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],"_",main,".csv"),row.names = F)
      if(file.exists(paste0(part,"docking/receptor_start/",df_all_systems$fin_name[i],".pdb"))){
        pdb<-read.pdb(paste0(part,"docking/receptor_start/",df_all_systems$fin_name[i],".pdb"))
        df_seq<-pdb$atom
        df_seq<-df_seq%>%filter(elety=="CA")
        df_seq<-df_seq%>%select(resno,resid, x, y, z)
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
        df_docking<-read.csv(paste0(part_start,"MD_analysis/docking/din/interaction_fin.csv"),stringsAsFactors = F)
        df_docking<-df_docking%>%mutate(system=as.character(system))
        df_docking<-df_docking%>%filter(system==df_all_systems$system_name[i])
        df_docking<-df_docking%>%select(resno, number_interactions, system, center, ligand, grops, grops_number, persent_interactions, aminoacids)
        df_seq<-left_join(df_seq,df_docking,by="resno")
        df_aligment<-read.csv(paste0(part_start,"fin_aligment.csv"),stringsAsFactors = F)
        df_seq<-left_join(df_seq,df_aligment,by="resno")
        df_seq<-df_seq%>%mutate(aminoacids=paste(resno,resid))
        
        write.csv(df_seq,paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],"_",main,".csv"),row.names = F)
        p_lenght<-ggplot(data = df_fin)+
          ggtitle(paste0("RMSD"))+
          labs(y = "RMSD (A)", x = "Time (ns)")+
          geom_line(aes(x = frame, y = s32,colour="303-322"))+
          geom_line(aes(x = frame, y = s38,colour="303-328"))+
          geom_line(aes(x = frame, y = s,colour="303-350"))+
          geom_line(aes(x = frame, y = s2,colour="322-350"))+
          geom_line(aes(x = frame, y = s8,colour="328-350"))+
          geom_line(aes(x = frame, y = sn,colour="322-328"))+
          theme_bw()+coord_flip()
        p_rmsd<-ggplot(data = df_fin)+
          ggtitle(paste0("RMSD"))+
          labs(y = "RMSD (A)", x = "Time (ns)")+
          geom_line(aes(x = frame, y = RMSD))+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$RMSD))
        p_sasa<-ggplot(data = df_fin)+
          labs(title=paste("SASA"),
               x = "Time(ns)", y = "SASA (A^2)") +
          geom_line(aes(x = frame,y=SASA))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$SASA))
  
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
          geom_hline(yintercept = median(df_fin$ramachadran))
  
        p_Total<-ggplot(data = df_fin)+
          labs(title=paste("Total energy"),
               x = "Time(ns)", y = "Total energy (kcal/mol)") +
          geom_line(aes(x = frame,y=Total))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$Total))
        p_VdW<-ggplot(data = df_fin)+
          labs(title=paste("VdW"),
               x = "Time(ns)", y = "VdW (kcal/mol)") +
          geom_line(aes(x = frame,y=VdW))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$VdW))
        p_rmsd_histo<-ggplot(data = df_fin)+
          ggtitle(paste0("RMSD"))+
          labs(x = "RMSD (A)")+
          geom_freqpoly(aes(x = RMSD))+
          theme_bw()+scale_fill_grey()+
          geom_text(x=median(df_RMSD$RMSD), y=100,label=round(median(df_RMSD$RMSD),digits = 1))+
          geom_vline(xintercept = median(df_RMSD$RMSD))
        p_rmsd_histo_EMD4<-ggplot(data = df_fin)+
          ggtitle(paste0("RMSD"))+
          labs(x = "RMSD (A)")+
          geom_freqpoly(aes(x = RMSD))+
          theme_bw()+scale_fill_grey()+
          geom_text(x=median(df_fin$RMSD), y=100,label=round(median(df_fin$RMSD),digits = 1))+
          geom_vline(xintercept = median(df_fin$RMSD))
        p_sasa_histo<-ggplot(data = df_fin)+
          labs(title="SASA", x = "SASA (A^2)") +
          geom_freqpoly(aes(x = SASA))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+
          geom_text(x=median(df_fin$SASA), y=100,label=round(median(df_fin$SASA),digits = 1))+
          geom_vline(xintercept = median(df_fin$SASA))
        p_Elec_histo<-ggplot(data = df_fin)+
          labs(title=paste("Electrostatic energy"),
               y = "", x = "Electrostatic energy (kcal/mol)") +
          geom_freqpoly(aes(x = Elec))+
          theme_bw()+
          geom_text(x=median(df_fin$Elec), y=100,label=round(median(df_fin$Elec),digits = 1))+
          geom_vline(xintercept = median(df_fin$Elec))
        p_ramachadran_histo <- ggplot(data = df_fin)+
          labs(title=paste("Ramachadran"),
               y = "", x = "Ramachadran") +
          geom_freqpoly(aes(x = ramachadran),bins=(max(df_fin$ramachadran)-min(df_fin$ramachadran)))+
          theme_bw()+#coord_flip()+
          geom_vline(xintercept = median(df_fin$ramachadran))
        p_Total_histo <- ggplot(data = df_fin)+
          labs(title=paste("Total energy"), y = "", x = "Total energy (kcal/mol)") +
          geom_freqpoly(aes(x = Total))+
          theme_bw()+
          geom_vline(xintercept = median(df_fin$Total))
        p_VdW_histo<-ggplot(data = df_fin)+
          labs(title=paste("VdW"), x = "VdW (kcal/mol)") +
          geom_freqpoly(aes(x =VdW))+
          theme_bw()+
          geom_text(x=median(df_fin$VdW), y=100,label=round(median(df_fin$VdW),digits = 1))+
          geom_vline(xintercept = median(df_fin$VdW))
        p_all<-plot_grid(p_sasa,       p_rmsd,      p_ramachadran,      p_Total,   
                         p_sasa_histo, p_rmsd_histo,p_ramachadran_histo,p_Total_histo,  
                         nrow=2,
                         rel_heights = c(4.5,1),
                         labels = c("A","B","C","D",
                                    "", "", "", "E"),align="hv")
        ggsave(p_all,   filename = paste0(part,"fin_plots/frame_plots/",df_all_systems$system_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 

        df_seqv<-df_seq%>%mutate(group_hbonds="NO")
        df_seqv$ligand[is.na(df_seqv$ligand)]<-"_NO"
        df_seqv$persent_interactions[is.na(df_seqv$persent_interactions)]<-0
        df_seqv<-df_seqv%>%mutate(group_docking="NO")
        df_seqv$grops_number[is.na(df_seqv$grops_number)]<-0
        df_seqv<-df_seqv%>%filter(grops_number<4)
        df_seqv$group_hbonds[df_seqv$hbonds>70]<-"water"
        df_seqv$group_docking[df_seqv$persent_interactions>70]<-"ligand"
#        df_seqv$group_docking[df_seqv$persent_interactions>50]<-df_seqv$ligand[df_seqv$persent_interactions>50]
        df_seqv<-df_seqv%>%mutate(group_interaction=paste(group_hbonds,group_docking))
        df_seqv<-df_seqv%>%filter(group_interaction!="NO NO")


        p_conserv<-ggplot(data = df_seqv)+
          geom_point(aes(x = resno, y = z,colour=persent_NaPi2b))+
          theme_bw()#+coord_flip()
        p_ramachadran<-ggplot(data = df_seqv)+
          geom_point(aes(x = resno, y = z,colour=ramachadran))+
          theme_bw()#+coord_flip()
        p_RMSF<-ggplot(data = df_seqv)+
          geom_point(aes(x = resno, y = z,colour=RMSF))+
          theme_bw()
        p_hbonds<-ggplot(data = df_seqv)+
          geom_point(aes(x = resno, y = z,colour=hbonds))+
          theme_bw()
        p_all<-plot_grid(p_conserv,p_ramachadran,       p_hbonds,      p_RMSF,     
                         nrow=4,  labels = c("A","B","C","D"),align="hv")
        ggsave(p_all,   filename = paste0(part,"fin_plots/str_plots/",df_all_systems$system_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
#        df_seqv<-df_seqv%>%mutate(is.na(ligand))  
        
        
        p_docking_Napi2b<-ggplot(data = df_seqv)+
          geom_vline(xintercept=c(157,158,159,160,431,432,433, 434))+
          geom_point(aes(x = resno, y = z,colour=persent_NaPi2b))+
          facet_grid(grops_number~ligand)+
          theme_bw()
        ggsave(p_docking_Napi2b,   filename = paste0(part,"fin_plots/docking_plots/NaPi2b_",df_all_systems$system_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
        p_docking_SLC34<-ggplot(data = df_seqv)+
          geom_vline(xintercept=c(157,158,159,160,431,432,433, 434))+
          geom_point(aes(x = resno, y = z,colour=persent_SLC34))+
          facet_grid(grops_number~ligand)+
          theme_bw()
        ggsave(p_docking_SLC34,   filename = paste0(part,"fin_plots/docking_plots/SLC34_",df_all_systems$system_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
        p_docking_SLC34<-ggplot(data = df_seqv)+
          geom_vline(xintercept=c(157,158,159,160,431,432,433, 434))+
          geom_point(aes(x = resno, y = z,colour=persent_SLC34))+
          facet_grid(group_interaction~ligand)+
          theme_bw()
        ggsave(p_docking_SLC34,   filename = paste0(part,"fin_plots/test_docking_plots/SLC34_",df_all_systems$system_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 30, units = c("cm"), dpi = 200 ) 
        p_docking_SLC34<-ggplot(data = df_seqv)+
          geom_vline(xintercept=c(157,158,159,160,431,432,433, 434))+
          geom_point(aes(x = resno, y = z,colour=persent_NaPi2b))+
          facet_grid(group_interaction~ligand)+
          theme_bw()
        ggsave(p_docking_SLC34,   filename = paste0(part,"fin_plots/test_docking_plots/NaPi2b_",df_all_systems$system_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",main,".png"), width = 60, height = 30, units = c("cm"), dpi = 200 ) 
       
 #       df_seqv<-df_seqv%>%select(aminoacids,center,ligand, group_interaction, persent_SLC34,persent_NaPi2b)
        
        df_seqv$persent_SLC34[df_seqv$persent_SLC34>80]<-"hight"
        df_seqv$persent_SLC34[df_seqv$persent_SLC34!="hight"]<-"low"
        df_seqv$persent_NaPi2b[df_seqv$persent_NaPi2b>80]<-"hight"
        df_seqv$persent_NaPi2b[df_seqv$persent_NaPi2b!="hight"]<-"low"
        
        df_seqv<-df_seqv%>%mutate(persent_SLC34=paste(center,ligand, group_interaction,persent_SLC34))
        df_seqv<-df_seqv%>%mutate(persent_NaPi2b=paste(center,ligand, group_interaction,persent_NaPi2b))
        df_seqv<-df_seqv%>%select(aminoacids,persent_SLC34,persent_NaPi2b)
        df_seqv<-unique(df_seqv)
        df_seqv_S<-df_seqv%>%group_by(persent_SLC34)%>%mutate(all=paste0(aminoacids,collapse = ", "))
        df_seqv_N<-df_seqv%>%group_by(persent_NaPi2b)%>%mutate(all=paste0(aminoacids,collapse = ", "))
        df_seqv_S<-ungroup(df_seqv_S)
        df_seqv_N<-ungroup(df_seqv_N)
        df_seqv_S<-df_seqv_S%>%select(persent_SLC34,all)
        df_seqv_N<-df_seqv_N%>%select(persent_NaPi2b,all)
        df_seqv_S<-unique(df_seqv_S)
        df_seqv_N<-unique(df_seqv_N)
        write.csv(df_seqv_S,paste0(part,"fin_data/str_data_sep/SLC34_",df_all_systems$fin_name[i],"_",main,".csv"),row.names = F)  
        write.csv(df_seqv_N,paste0(part,"fin_data/str_data_sep/NaPi2b_",df_all_systems$fin_name[i],"_",main,".csv"),row.names = F)
      }
    }
  }
}
