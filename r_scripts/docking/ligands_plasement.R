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
df_structure_surf<-read.csv(paste0(part,"docking/docking_first/din/df_merge_structure_log.csv"),stringsAsFactors = F)
df_structure_surf<-df_structure_surf%>%select(name.x,receptor,ligand,RMSD)
df_structure_surf<-df_structure_surf%>%mutate(RMSD=NA)
df_structure_surf<-unique(df_structure_surf)
df_structure_center<-df_structure_center%>%select(name.x,receptor,ligand,RMSD)
df_structure_center<-df_structure_center%>%mutate(RMSD=NA)
df_structure_center<-unique(df_structure_center)
df_structure<-left_join(df_structure_surf,df_structure_center,by=c("receptor","ligand","RMSD"))
i<-1
for (i in 1:nrow(df_structure)) {
  pdb_1<-read.pdb(paste0(part_start,"MD_analysis/docking/docking_first/din/structure_merged/",df_structure$name.x.x[i]))
  pdb_2<-read.pdb(paste0(part_start,"MD_analysis/docking/docking_first/din/structure_merged_center/",df_structure$name.x.y[i]))
  df_structure$RMSD[i]<-rmsd(pdb_1,pdb_2)
}

write.csv(df_structure,paste0(part,"fin_data/docking_statistic/comparition_surf_center_ligands.csv"),row.names = F)
#paste0(part,"fin_plots/energy_ligand_receptor_interactions.png")
df_structure<-read.csv(paste0(part,"fin_data/docking_statistic/comparition_surf_center_ligands.csv"),stringsAsFactors =  F)
df_temp<-df_structure%>%filter(RMSD<2.5)
v_surf<-df_structure$name.x.x[!(df_structure$name.x.x%in%df_temp$name.x.x)]
#print(length(v_surf))
v_center<-df_structure$name.x.x[(df_structure$name.x.x%in%df_temp$name.x.x)]
#print(length(v_center))

df_plot<-df_structure%>%filter(RMSD<10)
p<-ggplot(data=df_plot)+
  labs(x=expression(paste("RMSD ",
                          (ring(A)))))+
  geom_histogram(aes(x=RMSD,fill=NULL),binwidth=0.5)+
  theme_bw()
ggsave(p,filename = paste0(part,"fin_plots/docking_statistic/calibration_RMSD_ligand_allostaric_active_center.png"), width = 16, height = 12, units = c("cm"), dpi = 1000 ) 




df_structure_center<-read.csv(paste0(part,"docking/docking_first/din/df_merge_structure_log_center.csv"),stringsAsFactors = F)
df_structure_center<-df_structure_center%>%select(name.x,receptor,ligand,affinity)
df_structure_center<-df_structure_center%>%mutate(folder_interactions="interaction_center")
df_structure_center<-df_structure_center%>%mutate(folder_structure="structure_merged_center")
df_structure_center<-df_structure_center%>%mutate(type="active center")
nrow(df_structure_center)
df_structure_surf<-read.csv(paste0(part,"docking/docking_first/din/df_merge_structure_log.csv"),stringsAsFactors = F)
df_structure_surf<-df_structure_surf%>%select(name.x,receptor,ligand,affinity)
df_structure_surf<-df_structure_surf%>%mutate(folder_interactions="interaction_serf")
df_structure_surf<-df_structure_surf%>%mutate(folder_structure="structure_merged")
df_structure_surf<-df_structure_surf%>%mutate(type="allosteric interactions")

#nrow(df_structure_surf)
df_structure_surf$type[df_structure_surf$name.x%in%v_center]<-"active center"
#unique(df_structure_surf$type)
df_ligands_type<-rbind(df_structure_surf,df_structure_center)

df_ligands_type<-df_ligands_type%>%mutate(parset=paste(receptor,ligand,type))
df_ligands_type<-left_join(df_ligands_type,df_all_systems,by=c("receptor"="fin_name"))
write.csv(df_ligands_type, paste0(part,"fin_data/docking_statistic/energy_receptor_ligand.csv"),row.names = F)

kruskal <- df_ligands_type %>% kruskal_test(affinity ~ parset)


pwc2 <- df_ligands_type %>% 
  dunn_test(affinity ~ parset, p.adjust.method = "bonferroni") 

df_convert<-df_ligands_type%>%select(receptor,ligand,type,parset)
df_convert<-unique(df_convert)
#colnames(df_ligands_type)

pwc2 <- pwc2 %>% add_xy_position(x = "parset")
pwc2<-left_join(pwc2,df_convert,by=c("group1"="parset"))
pwc2<-left_join(pwc2,df_convert,by=c("group2"="parset"))
pwc2<-pwc2%>%filter(p.adj<0.05)
pwc2<-pwc2%>%filter(receptor.x==receptor.y)%>%filter(ligand.x== ligand.y)
pwc2<-pwc2%>%mutate(receptor=receptor.x)
pwc2<-pwc2%>%mutate(ligand=ligand.x)
pwc2<-pwc2%>%select(group1,group2,  statistic,p,p.adj,p.adj.signif,receptor,ligand)
pwc2<-left_join(pwc2,df_all_systems,by=c("receptor"="fin_name"))


p<-ggplot(data=df_ligands_type)+
  labs(x="Affinity, kcal/mol")+
  geom_boxplot(aes(x=affinity,y=type))+
  geom_text(aes(x=1.5,y=1,5,label=format(p.adj, digits=3, big.mark   = ",")),data=pwc2)+
  facet_grid(system_name~ligand)+
  theme_bw()
ggsave(p,filename = paste0(part,"fin_plots/docking_statistic/type_energy_ligand_receptor_interactions.png"), width = 20, height = 16, units = c("cm"), dpi = 1000 ) 

