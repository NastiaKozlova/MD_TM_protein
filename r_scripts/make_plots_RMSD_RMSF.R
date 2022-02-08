part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(dplyr)
library(Peptides)
library(ggplot2)
library(cowplot)
library(bio3d)

test_10<-seq(from=0,to=1000,by=10)
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
if (!dir.exists(paste0(part,"fin_data/str_data"))) {dir.create(paste0(part,"fin_data/str_data"))}
if (!dir.exists(paste0(part,"fin_data/aminoacids_interactions/"))){dir.create(paste0(part,"fin_data/aminoacids_interactions/"))}

if (!dir.exists(paste0(part,"fin_plots"))) {dir.create(paste0(part,"fin_plots"))}
if (!dir.exists(paste0(part,"fin_plots/str_plots"))) {dir.create(paste0(part,"fin_plots/str_plots"))}
if (!dir.exists(paste0(part,"fin_plots/frame_plots"))) {dir.create(paste0(part,"fin_plots/frame_plots"))}
if (!dir.exists(paste0(part,"fin_plots/docking_plots"))) {dir.create(paste0(part,"fin_plots/docking_plots"))}

i<-1
main<-main_part[1]
for (i in 1:nrow(df_all_systems)) {
  for (main in main_part) {
    if (file.exists(paste0(parta,df_all_systems$fin_name[i],"/","din/Energy/protein_",main,".txt"))){
      if (file.exists(paste0(part,"din/",df_all_systems$fin_name[i],"/",main,"_time_Ramachadran.csv"))){
        df_all_systems$Progress[i]<-"done"
        df_second_structure<-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/Second_structure_",main,".csv"),stringsAsFactors = F)
        df_second_structure<-read.csv(paste0(part,"din/",df_all_systems$fin_name[i],"/Second_structure_",main,".csv"),stringsAsFactors = F)
        
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
        df_fin<-full_join(df_ramachadran,df_Energy_protein,by=c("frame"="Frame"))
        df_fin$time<-NULL
        df_fin$Time<-NULL
        df_fin<-left_join(df_fin,df_Energy_protein_water,by=c("frame"="Frame"))
        df_fin$Time<-NULL
        df_fin$time<-NULL
        df_fin<-left_join(df_fin,df_Energy_protein_lipids,by=c("frame"="Frame"))
        df_fin$Time<-NULL
        df_fin$time<-NULL
        df_SASA<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/SASA/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
         colnames(df_SASA)<-c("frame","SASA_protein")
        df_fin<-left_join(df_fin,df_SASA,by=c("frame"))
        df_fin$Time<-NULL
        df_fin$time<-NULL
        df_RMSD_protein<-read.table(paste0(parta,df_all_systems$fin_name[i],"/","din/RMSD/",main,".txt"), sep="", header=F, na.strings ="", stringsAsFactors= F)
        colnames(df_RMSD_protein)<-c("frame","RMSD_protein")
        df_fin<-left_join(df_fin,df_RMSD_protein,by=c("frame"))
        df_fin$Time<-NULL
        df_fin$time<-NULL 
        df_fin<-df_fin%>%mutate(frame=frame/10)
        write.csv(df_fin,paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],".csv"),row.names = F)
        p_second<-ggplot(data = df_second_structure)+
          ggtitle(paste0("Second structure ",main))+
          labs(x = "Number of aminoasids", y = "Time (ns)")+
          geom_rect(aes(xmin = start, ymin = level_min, xmax= finish, ymax = level_max, colour = type,fill=type))+
          scale_color_manual(values = v_pallete)+ scale_fill_manual(values = v_pallete)+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          scale_y_continuous(breaks = test_10, minor_breaks = NULL)+
          theme_bw()+theme(legend.position = "top")+guides(fill="none",colour="none")
        p_RMSD<-ggplot(data = df_fin)+
          ggtitle(paste0("RMSD"))+
          labs(y = "RMSD (A)", x = "Time (ns)")+
          geom_line(aes(x = frame, y = RMSD_protein))+
          theme_bw()+coord_flip()+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          geom_hline(yintercept = median(df_fin$RMSD_protein,na.rm = T))
        p_sasa<-ggplot(data = df_fin)+
          labs(title=paste("SASA"),
               x = "Time(ns)", y = "SASA (A^2)") +
          geom_line(aes(x = frame,y=SASA_protein))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$SASA_protein,na.rm = T))
        
        p_Elec<-ggplot(data = df_fin)+
          labs(title=paste("Electrostatic energy"),
               x = "Time(ns)", y = "Electrostatic energy (kcal/mol)") +
          geom_line(aes(x = frame,y=protein_Elec))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$protein_Elec))
        p_ramachadran<-ggplot(data = df_fin)+
          labs(title=paste("Ramachadran"), x = "Time(ns)", y = "Ramachadran") +
          geom_line(aes(x = frame,y=ramachadran))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$ramachadran,na.rm = T))
        
        p_Total<-ggplot(data = df_fin)+
          labs(title=paste("Total energy"),
               x = "Time(ns)", y = "Total energy (kcal/mol)") +
          geom_line(aes(x = frame,y=protein_Total))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$protein_Total,na.rm = T))
        p_VdW<-ggplot(data = df_fin)+
          labs(title=paste("VdW"),
               x = "Time(ns)", y = "VdW (kcal/mol)") +
          geom_line(aes(x = frame,y=protein_VdW))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+coord_flip()+
          geom_hline(yintercept = median(df_fin$protein_VdW,na.rm = T))
        p_rmsd_histo<-ggplot(data = df_fin)+
          ggtitle(paste0("RMSD"))+
          labs(x = "RMSD (A)")+
          geom_freqpoly(aes(x = RMSD_protein))+
          theme_bw()+scale_fill_grey()+
          geom_text(x=median(df_fin$RMSD_protein), y=20,label=round(median(df_fin$RMSD_protein,na.rm = T),digits = 1))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          geom_vline(xintercept = median(df_fin$RMSD_protein,na.rm = T))
        p_sasa_histo<-ggplot(data = df_fin)+
          labs(title="SASA", x = "SASA (A^2)") +
          geom_freqpoly(aes(x = SASA_protein))+
          theme_bw()+
          geom_text(x=median(df_fin$SASA_protein), y=20,label=round(median(df_fin$SASA_protein,na.rm = T),digits = 1))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          geom_vline(xintercept = median(df_fin$SASA_protein,na.rm = T))
        p_Elec_histo<-ggplot(data = df_fin)+
          labs(title=paste("Electrostatic energy"),
               y = "", x = "Electrostatic energy (kcal/mol)") +
          geom_freqpoly(aes(x = protein_Elec))+
          theme_bw()+
          geom_text(x=median(df_fin$protein_Elec,na.rm = T), y=20,label=round(median(df_fin$protein_Elec,na.rm = T),digits = 1))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          geom_vline(xintercept = median(df_fin$protein_Elec,na.rm = T))
        p_ramachadran_histo <- ggplot(data = df_fin)+
          labs(title=paste("Ramachadran"),
               y = "", x = "Ramachadran") +
          geom_freqpoly(aes(x = ramachadran),bins=(max(df_fin$ramachadran)-min(df_fin$ramachadran)))+
          theme_bw()+#coord_flip()+
          geom_text(x=median(df_fin$ramachadran,na.rm = T), y=20,label=round(median(df_fin$ramachadran,na.rm = T),digits = 1))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          geom_vline(xintercept = median(df_fin$ramachadran,na.rm = T))
        p_Total_histo <- ggplot(data = df_fin)+
          labs(title=paste("Total energy"), y = "", x = "Total energy (kcal/mol)") +
          geom_freqpoly(aes(x = protein_Total))+
          theme_bw()+
          geom_vline(xintercept = median(df_fin$protein_Total,na.rm = T))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          geom_text(x=median(df_fin$protein_Total,na.rm = T), y=20,label=round(median(df_fin$protein_Total,na.rm = T),digits = 1))
        p_VdW_histo<-ggplot(data = df_fin)+
          labs(title=paste("VdW"), x = "VdW (kcal/mol)") +
          geom_freqpoly(aes(x =protein_VdW))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()+
          geom_text(x=median(df_fin$protein_VdW,na.rm = T), y=20,label=round(median(df_fin$protein_VdW,na.rm = T),digits = 1))+
          geom_vline(xintercept = median(df_fin$protein_VdW,na.rm = T))
        p_rmsf<-ggplot(data = df_RMSF)+
          ggtitle(paste0("RMSF"))+
          labs(x = "Number of aminoasids", y = "RMSF (A)")+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          #labs(x = "Номер аминокислоты", y = "RMSF (A)")+
          geom_line(aes(x = Resid, y = RMSF))+
          theme_bw()
        p_all<-plot_grid(p_sasa,       p_RMSD,      p_ramachadran,      p_Total,   p_second,
                         p_sasa_histo, p_rmsd_histo,p_ramachadran_histo,p_Total_histo,  p_rmsf,
                         ncol=5,
                         rel_heights = c(4.5,1),rel_widths = c(1,1,1,1,5),
                         labels = c("A","B","C","D","E",
                                    "", "", "", "",""),align="hv")
        ggsave(p_all,   filename = paste0(part,"fin_plots/frame_plots/",df_all_systems$fin_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],"_",df_all_systems$Structure[i],".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
        
        
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
        p_conserv<-ggplot(data = df_seq)+
          geom_point(aes(x = resno, y = conservative,colour=topology))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()#+coord_flip()
        p_ramachadran<-ggplot(data = df_seq)+
          geom_point(aes(x = resno, y = ramachadran,colour=topology))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()#+coord_flip()
        p_RMSF<-ggplot(data = df_seq)+
          geom_point(aes(x = resno, y = RMSF,colour=topology))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()
        p_hbonds<-ggplot(data = df_seq)+
          geom_point(aes(x = resno, y = hbonds,colour=topology))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()
        p_importance_interprotein_intractions<-ggplot(data = df_seq)+
          geom_point(aes(x = resno, y = importance_interprotein_intractions,colour=topology))+
          scale_x_continuous(breaks = test_10, labels =  test_10)+
          theme_bw()
        p_all<-plot_grid(p_conserv,p_ramachadran,       p_hbonds,      p_RMSF,     p_importance_interprotein_intractions,
                         nrow=5,  labels = c("A","B","C","D","E"),align="hv")
        ggsave(p_all,   filename = paste0(part,"fin_plots/str_plots/Membrane_",df_all_systems$Membrane[i],
                                          "_pH_",df_all_systems$fin_name[i],"_Structure_",df_all_systems$Structure[i],
                                          "_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],
                                          "_",df_all_systems$fin_name[i],".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
      }
    }
  }
}
if (!dir.exists(paste0(part,"fin_data/docking_data"))) {dir.create(paste0(part,"fin_data/docking_data"))}
i<-1 
j<-1
for (i in 1:nrow(df_all_systems)) {
  df_str<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors = F)
  df_system<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/df_all.csv"),stringsAsFactors =  F)
  df_system<-df_system%>%select(receptor,ligand)
  df_system<-unique(df_system)
  df_system<-df_system%>%filter(receptor==df_all_systems$fin_name[i])
  df_system<-left_join(df_system,df_all_systems,by=c("receptor"="fin_name"))
  df_docking<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/din/interaction_fin/",
                            df_system$receptor[1],"_",df_system$ligand[1],".csv"),stringsAsFactors = F)
  if(nrow(df_system)>1){
  for (j in 2:nrow(df_system)) {
    df_docking_add<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/din/interaction_fin/",
                                    df_system$receptor[j],"_",df_system$ligand[j],".csv"),stringsAsFactors = F)
    df_docking<-rbind(df_docking,df_docking_add)
    df_docking<-unique(df_docking)
    }
  }

  df_docking<-df_docking%>%mutate(persent_receptor_ligand_interaction=total_persent_interactions)
  df_docking<-df_docking%>%select(resid,resno,ligand,size_of_group,persent_receptor_ligand_interaction)
  df_docking<-unique(df_docking)
  df_str<-left_join(df_str,df_docking,by=c("resno","resid"))
  write.csv(df_str,paste0(part,"fin_data/docking_data/",df_all_systems$fin_name[i],".csv"),row.names = F )
  df_str<-df_str%>%filter(persent_receptor_ligand_interaction>0)
  p_docking<-ggplot(data = df_str)+#geom_rect(aes(xmin=seq_beg, xmax=seq_end, ymin=-Inf, ymax=Inf, fill=topology),data=df_topology)+
    geom_text(aes(x = resno, y = z,colour=conservative,label=persent_receptor_ligand_interaction))+
    theme_bw()+ facet_grid(size_of_group~ligand)
  ggsave(p_docking,   filename = paste0(part,"fin_plots/docking_plots/",df_all_systems$fin_name[i],"_disulfid_bonds_",df_all_systems$disulfid_bonds[i],"_glyco_",df_all_systems$Glyco[i],"_",df_all_systems$Membrane[i],".png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 
  
}


df_seq<-read.csv(paste0(part,"fin_data/docking_data/",df_all_systems$fin_name[1],".csv"),stringsAsFactors =  F)
df_seq<-df_seq%>%mutate(receptor=df_all_systems$fin_name[1])
print(nrow(df_seq))
if (nrow(df_all_systems)>1) {
  for (i in 2:nrow(df_all_systems)) {
    df_seq_add<-read.csv(paste0(part,"fin_data/docking_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
    df_seq_add<-df_seq_add%>%mutate(receptor=df_all_systems$fin_name[i])
    print(nrow(df_seq_add))
    df_seq<-rbind(df_seq,df_seq_add)
  }
}
print(nrow(df_seq))
df_seq_test<-df_seq%>%group_by(receptor)%>%mutate(size=n())
df_seq_test<-df_seq_test%>%select(receptor,size)
df_seq_test<-unique(df_seq_test)
df_seq<-df_seq%>%filter(persent_receptor_ligand_interaction>0)
df_seq_test1<-df_seq%>%group_by(receptor)%>%mutate(size=n())
df_seq_test1<-df_seq_test1%>%select(receptor,size)
df_seq_test1<-unique(df_seq_test1)
#df_seq<-df_seq%>%filter(ramachadran==0)
#df_seq<-df_seq%>%filter(conservative>75)
#df_seq<-df_seq%>%filter(hbonds>50)
p_docking<-ggplot(data = df_seq)+#geom_rect(aes(xmin=seq_beg, xmax=seq_end, ymin=-Inf, ymax=Inf, fill=type),data=df_topology)+
  geom_text(aes(x = resno, y = resid,colour=conservative,label=amino))+
  theme_bw()+ facet_grid(ligand~receptor)
ggsave(p_docking,   filename = paste0(part,"fin_plots/docking_fin_intercations.png"), width = 60, height = 40, units = c("cm"), dpi = 200 ) 

df_all_systems<-df_all_systems%>%mutate(system_name=paste(Membrane,Structure))
df_structure_RMSD<-read.csv(paste0(part,"docking/docking_first/din/df_merge_structure_log.csv"),stringsAsFactors = F)
df_structure_RMSD<-left_join(df_structure_RMSD,df_all_systems,by=c("receptor"="fin_name"))

v_min<-quantile(df_structure_RMSD$affinity,0.025)
v_max<-quantile(df_structure_RMSD$affinity,0.975)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%filter(affinity>v_min)
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%filter(affinity<v_max)


a<-seq(from=min(df_structure_RMSD$affinity),to=max(df_structure_RMSD$affinity),by=0.1)
df_structure_RMSD<-df_structure_RMSD%>%mutate(number=as.character(number))
p<-ggplot(data=df_structure_RMSD)+geom_freqpoly(aes(x=affinity,colour=number),binwidth=0.1)+facet_grid(system_name~ligand)+
  scale_x_continuous(breaks=a,labels=a)+theme_bw()+guides(color = "none", size = "none")+
  labs(x="Affinity, kcal/mol")

ggsave(p,filename = paste0(part,"fin_plots/energy_ligand_receptor_interactions.png"), width = 12, height = 12, units = c("cm"), dpi = 1000 ) 


