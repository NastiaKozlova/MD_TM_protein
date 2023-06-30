part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(dplyr)
library(Peptides)
library(ggplot2)
library(cowplot)
library(bio3d)
library(rstatix)


test_10<-seq(from=0,to=1000,by=10)
v_pallete<-c("sheet"="#BBBBBB","helix"="#333333")
df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors =  F)
part<-paste0(part_start,"MD_analysis/")
if (!dir.exists(paste0(part,"fin_data/docking_statistic"))) {dir.create(paste0(part,"fin_data/docking_statistic"))}
#parta<-paste0(part_start,"MD/")
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)


df_all_systems$system_name<-as.character(df_all_systems$system_name)
#df_all_systems<-df_all_systems%>%mutate(Progress=="notdone")
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
#\    df_all_systems$Membrane[i],"_",df_all_systems$Structure[i]
df_all_systems<-df_all_systems%>%mutate(system_name=paste(Membrane,Structure))
v_part<-list.files(paste0(part,"din"))
df_structure_center<-read.csv(paste0(part,"docking/docking_first/din/df_merge_structure_log_center.csv"),stringsAsFactors = F)
#df_structure_center<-df_structure_center%>%select(name.x,receptor,ligand,affinity)
df_structure_center<-df_structure_center%>%filter(models.x==models.y)
df_structure_center<-df_structure_center%>%filter(name.x==name.y)
df_structure_center<-df_structure_center%>%select(name.x,receptor,ligand,affinity)
df_structure_center<-unique(df_structure_center)
#unique(df_structure_center$name.x)
df_structure_center<-df_structure_center%>%mutate(type="center")

df_structure_center<-unique(df_structure_center)

df_structure_center<-df_structure_center%>%mutate(x=NA)
df_structure_center<-df_structure_center%>%mutate(y=NA)
df_structure_center<-df_structure_center%>%mutate(z=NA)
i<-2
#length(unique(df_structure_center$name.x))
for (i in 1:nrow(df_structure_center)) {
  pdb<-read.pdb(paste0(part,"docking/docking_first/din/str_fin/",df_structure_center$name.x[i]))
  df_pdb<-pdb$atom
  df_structure_center$x[i]<-median(df_pdb$x)
  df_structure_center$y[i]<-median(df_pdb$y)
  df_structure_center$z[i]<-median(df_pdb$z)
  
  df_structure_center$z_min[i]<-min(df_pdb$z)
  df_structure_center$z_max[i]<-max(df_pdb$z)
}


df_structure_surf<-read.csv(paste0(part,"docking/docking_first/din/df_merge_structure_log.csv"),stringsAsFactors = F)
df_structure_surf<-df_structure_surf%>%filter(models.x==models.y)
df_structure_surf<-df_structure_surf%>%mutate(name.y=name.x)

#df_structure_surf<-df_structure_surf%>%filter(RMSD==0)
df_structure_surf<-df_structure_surf%>%select(name.x,receptor,ligand,affinity)

df_structure_surf<-unique(df_structure_surf)

df_structure_surf<-df_structure_surf%>%mutate(type="surf")



df_structure_surf<-df_structure_surf%>%mutate(x=NA)
df_structure_surf<-df_structure_surf%>%mutate(y=NA)
df_structure_surf<-df_structure_surf%>%mutate(z=NA)
df_structure_surf<-df_structure_surf%>%mutate(z_min=NA)
df_structure_surf<-df_structure_surf%>%mutate(z_max=NA)

for (i in 1:nrow(df_structure_surf)) {
  pdb<-read.pdb(paste0(part,"docking/docking_first/din/str_fin/",df_structure_surf$name.x[i]))
  df_pdb<-pdb$atom
  df_structure_surf$x[i]<-median(df_pdb$x)
  df_structure_surf$y[i]<-median(df_pdb$y)
  df_structure_surf$z[i]<-median(df_pdb$z)
  df_structure_surf$z_min[i]<-min(df_pdb$z)
  df_structure_surf$z_max[i]<-max(df_pdb$z)
}
df_structure<-rbind(df_structure_center,df_structure_surf)


df_structure<-left_join(df_structure,df_all_systems,by=c("receptor"="fin_name"))
df_structure<-df_structure%>%group_by(receptor,ligand)%>%arrange(desc(z))%>%mutate(structure_order=c(1:n()))
df_structure<-df_structure[df_structure$system_name%in%c("POPE WT","POPG Inverted"),]
#df_structure<-df_structure[df_structure$system_name%in%c("POPE WT"),]
df_structure<-df_structure%>%filter(ligand!="D-Glucose")
p<-ggplot(data=df_structure)+
  geom_segment(aes(x=z_min,xend=z_max,y = affinity,yend=affinity,colour=type))+
#  geom_text(aes(x=z,y=affinity,colour=type,label=structure_order))+
#  geom_line(aes(x=structure_order,y=z))+
#  geom_point(aes(x=structure_order,y=z))+
  
  geom_smooth(aes(x=z,y=affinity))+
  facet_grid(ligand~system_name)+
  
  theme_bw()+
  geom_hline(yintercept = 0)

ggsave(p,filename = paste0(part,"fin_plots/docking_statistic/interaction_between_energy_z_lingand_receptor_all.png"), width = 24, height = 20, units = c("cm"), dpi = 200 )
df_structure<-df_structure%>%group_by(receptor,ligand)%>%arrange(desc(z))%>%mutate(structure_order=c(1:n()))
write.csv(df_structure,paste0(part,"docking/docking_first/din/df_structure.csv"),row.names = F)
