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
i<-1


df_frame_data<-read.csv(paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
df_seq_data<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
df_ring<-read.csv(paste0(part,"fin_data/aminoacids_interactions/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)

df_frame_data<-df_frame_data%>%mutate(fin_name=df_all_systems$fin_name[1])
df_seq_data<-df_seq_data%>%mutate(fin_name=df_all_systems$fin_name[1])
df_ring<-df_ring%>%mutate(fin_name=df_all_systems$fin_name[1])
for (i in 2:nrow(df_all_systems)) {
    df_frame_data_add<-read.csv(paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
    df_seq_data_add<-read.csv(paste0(part,"fin_data/str_data/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
    df_ring_add<-read.csv(paste0(part,"fin_data/aminoacids_interactions/",df_all_systems$fin_name[i],".csv"),stringsAsFactors =  F)
    
    df_frame_data_add<-df_frame_data_add%>%mutate(fin_name=df_all_systems$fin_name[i])
    df_seq_data_add<-df_seq_data_add%>%mutate(fin_name=df_all_systems$fin_name[i])
    df_ring_add<-df_ring_add%>%mutate(fin_name=df_all_systems$fin_name[i])
    
    df_frame_data<-rbind(df_frame_data,df_frame_data_add)
    df_seq_data<-rbind(df_seq_data,df_seq_data_add)
    df_ring<-rbind(df_ring,df_ring_add)
}
df_frame_data<-left_join(df_frame_data,df_all_systems,by = c("fin_name"))
#df_frame_data<-df_frame_data%>%mutate(system_name=paste(Membrane,Structure))
df_summary<-df_frame_data%>%group_by(fin_name)%>%summarise(max_frame=max(frame))
df_frame_data<-df_frame_data%>%filter(frame<=min(df_summary$max_frame))
#df_frame_data<-left_join(df_frame_data,df_all_systems,by="fin_name")

p_RMSD<-ggplot(data = df_frame_data)+
  ggtitle(paste0("RMSD"))+
  labs(y = "RMSD (A)", x = "Time (ns)")+
  geom_line(aes(x = frame, y = RMSD_protein))+
  facet_grid(Structure~Membrane)+#,colour=Structure,linetype=Membrane))+
  theme_bw()+
  theme(legend.position = "bottom")+
 #coord_flip()+
  scale_x_continuous(breaks = test_10, labels =  test_10)#+

ggsave(p_RMSD,   filename = paste0(part,"fin_plots/frame_statisitc/RMSD_graph.png"), width = 20, height = 15, units = c("cm"), dpi = 600 ) 

df_frame_data<-df_frame_data%>%mutate(protein_protein_Total=protein_Total-protein_water_Total-protein_lipid_Total)
#df_frame_data<-left_join(df_frame_data,df_all_systems,by=c("system"="fin_name"))
#df_frame_data<-df_frame_data%>%filter(frame>10)
df_frame_data<-df_frame_data%>%mutate(system=paste0(Structure,"_",Membrane))
p_protein_lipid_histo<-ggplot(data = df_frame_data)+
  labs(x = "Energy (kcal/mol)",title="Energy of protein lipid \ninteractions",) +theme_bw()+
  geom_freqpoly(aes(x =protein_lipid_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
p_protein_water_histo<-ggplot(data = df_frame_data)+
  labs(x = "Energy (kcal/mol)",title="Energy of protein water\n interactions") +theme_bw()+
  geom_freqpoly(aes(x =protein_water_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
p_protein_protein_histo<-ggplot(data = df_frame_data)+
  labs(x = "Energy (kcal/mol)",title="Energy of \ninterprotein interactions") +theme_bw()+
  geom_freqpoly(aes(x =protein_protein_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
p_protein_histo<-ggplot(data = df_frame_data)+
  labs(x = "Energy (kcal/mol)",title="Total \nprotein energy") +theme_bw()+
  geom_freqpoly(aes(x =protein_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
legend_test<-get_legend(p_protein_histo)

p_all<-plot_grid(p_protein_lipid_histo+theme(legend.position = "none"),
                 p_protein_water_histo+theme(legend.position = "none"),  
                 p_protein_protein_histo+theme(legend.position = "none"),
                 p_protein_histo+theme(legend.position = "none"),
                 rel_heights=c(1,1),
                 nrow=2,  labels = c("A","B","C","D"),align="hv",ncol = 2,axis="bt")
p_test<-plot_grid(p_all,legend_test,nrow=2,rel_heights=c(1,0.1))
ggsave(p_test,filename = paste0(part,"fin_plots/frame_statisitc/protein_energy_compation.png"), width = 15, height = 12, units = c("cm"), dpi = 1000 ) 
