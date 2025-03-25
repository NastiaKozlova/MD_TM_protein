part_start <- commandArgs(trailingOnly=TRUE)
library(dplyr)
library(ggplot2)
library(bio3d)
seq_name<-list.files(paste0(part_start,"start/sequence/"))
seqv<-read.fasta(paste0(part_start,"start/sequence/",seq_name))
#df_systems<-read.csv("start/all_systems.csv",stringsAsFactors = F)
#df_systems<-df_systems%>%mutate(system_name=paste(Membrane,Structure))
v_pallete<-c("CD"="#ffc785","ECD"="#A5F6B1","TMD"="#f6a5c6")
v_seq<-length(seqv$ali)
v_sep<-seq(from=0,to=v_seq,by=10)
part_name<-list.files(paste0(part_start,"MD_analysis/"))
df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors = F)
df_interaction<-read.csv(paste0(part_start,
                                "MD_analysis/fin_data/docking_statistic/path_aminoacid_interactions.csv"),stringsAsFactors = F)

df_interaction<-df_interaction%>%filter(persent_interactions>90)
df_interaction<-df_interaction%>%filter(affinity<0)
#df_interaction<-left_join(df_interaction,df_systems,by="system_name")
df_interaction<-df_interaction%>%mutate(z.y=orientation*z.y)
#df_interaction<-df_interaction%>%mutate(z.x=orientation*z.x)
p<-ggplot(data=df_interaction)+
  geom_rect(aes(xmin = seq_beg-0.5, ymin = -Inf,
                xmax = seq_end+0.5, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
  geom_rect(aes(xmin = resno-0.5,xmax=resno+0.5,ymin=z.y-.1,ymax=z.y+.1,colour=type.y))+#,fill=type.y))+
  geom_hline(yintercept =c(-20,0,20) )+
  facet_grid(system_name~ligand)+
  theme_bw()
ggsave(p,filename = paste0(part_start,
                           "MD_analysis/fin_plots/docking_statistic/interaction_between_interactions_z_lingand_placement.png"),
       width = 30, height = 20, units = c("cm"), dpi = 200 )
p<-ggplot(data=df_interaction)+
  geom_rect(aes(xmin = seq_beg-0.5, ymin = -Inf,
                xmax = seq_end+0.5, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
  geom_rect(aes(xmin = resno-0.5,xmax=resno+0.5,ymin=z.y-.1,ymax=z.y+.1,colour=affinity))+#,fill=type.y))+
  geom_hline(yintercept =c(-20,0,20) )+
  facet_grid(system_name~ligand)+
  theme_bw()
ggsave(p,filename = paste0(part_start,
                           "MD_analysis/fin_plots/docking_statistic/interaction_between_interactions_z_receptor_placement.png"),
       width = 30, height = 20, units = c("cm"), dpi = 200 )

v_system_name<-unique(df_interaction$system_name)
for (i in 1:length(v_system_name)) {
  df_interaction_sep<-df_interaction%>%filter(system_name==v_system_name[i])
  p<-ggplot(data=df_interaction_sep)+
    
    geom_rect(aes(xmin = seq_beg-0.5, ymin = -Inf,
                  xmax = seq_end+0.5, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
    geom_rect(aes(xmin = resno-0.5,xmax=resno+0.5,ymin=z.y-.1,ymax=z.y+.1,colour=type.y))+#,fill=type.y))+
    scale_x_continuous(breaks = v_sep,labels =v_sep )+ 
    scale_fill_manual(values=v_pallete)+
    geom_hline(yintercept =c(-20,0,20) )+
    facet_grid(ligand~.)+
    theme_bw()
  ggsave(p,filename = paste0(part_start,
                             "MD_analysis/fin_plots/docking_statistic/interaction_between_interactions_",
                             v_system_name[i],"_z_lingand_placement.png"),
         width = 30, height = 20, units = c("cm"), dpi = 200 )
  p<-ggplot(data=df_interaction_sep)+
    geom_rect(aes(xmin = seq_beg-0.5, ymin = -Inf,
                  xmax = seq_end+0.5, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
    geom_rect(aes(xmin = resno-0.5,xmax=resno+0.5,ymin=z.y-.1,ymax=z.y+.1,colour=affinity))+#,fill=type.y))+
    scale_x_continuous(breaks = v_sep,labels =v_sep )+ 
    scale_fill_manual(values=v_pallete)+
    geom_hline(yintercept =c(-20,0,20) )+
    facet_grid(ligand~.)+
    theme_bw()
  ggsave(p,filename = paste0(part_start,
                             "MD_analysis/fin_plots/docking_statistic/interaction_between_interactions_",
                             v_system_name[i],"_z_receptor_placement.png"),
         width = 30, height = 20, units = c("cm"), dpi = 200 )
}
df_interaction_test<-df_interaction#[df_interaction$resno%in%c(240,243,246,247,250,254),]
v_sep<-seq(from=0,to=v_seq,by=2)
df_interaction_test<-df_interaction_test[df_interaction_test$system_name%in%c("POPE WT","POPG Inverted"),]
df_interaction_test<-df_interaction_test[df_interaction_test$ligand%in%c("POPE","POPG"),]
p<-ggplot(data=df_interaction_test)+
  geom_rect(aes(xmin = seq_beg-0.5, ymin = -Inf,
                xmax = seq_end+0.5, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
  geom_rect(aes(xmin = resno-0.5,xmax=resno+0.5,ymin=z.y-.1,ymax=z.y+.1,colour=affinity))+#,fill=type.y))+
  scale_x_continuous(breaks = v_sep,labels =v_sep,limits=c(235,260) )+ 
  scale_fill_manual(values=v_pallete)+
  
  geom_hline(yintercept =c(-20,0,20) )+
  facet_grid(system_name~ligand)+
  theme_bw()
#ggsave(p,filename = paste0(part_start,
#                           "MD_analysis/fin_plots/docking_statistic/interaction_between_interactions_",
#                           v_system_name[i],"_z_receptor_placement.png"),
#       width = 30, height = 20, units = c("cm"), dpi = 200 )
ggsave(p,filename = paste0(part_start,
                           "MD_analysis/fin_plots/docking_statistic/tost.png"),
       width = 30, height = 20, units = c("cm"), dpi = 200 )
