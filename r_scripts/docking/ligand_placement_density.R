part_analysis <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(ggplot2)
library(dplyr)
v_z<-seq(from=-60,to=60,by=1)
df_pore<-data.frame(matrix(ncol = 2,nrow = length(v_z)))
colnames(df_pore)<-c("z","size")
df_pore$z<-v_z
df_pore$size<-0

part<-paste0(part_analysis,"din/")
setwd(part)
df_structure_RMSD<-read.csv("df_merge_structure_log.csv",stringsAsFactors = F)

#df_structure_RMSD<-df_structure_RMSD%>%filter(RMSD==0)
df_structure_RMSD<-df_structure_RMSD%>%select(name.x,receptor,ligand)
df_structure_RMSD<-unique(df_structure_RMSD)
if(dir.exists(paste0(part,"ligand_placement_density"))) {system(command = paste0("rm -r ",part,"ligand_placement_density"),ignore.stdout=T,wait = T)}
if (!dir.exists("ligand_placement_density")) {dir.create("ligand_placement_density")}
if (!dir.exists("ligand_placement_density_merge")) {dir.create("ligand_placement_density_merge")}
if (!dir.exists("ligand_placement_density_merge_data_plot")) {dir.create("ligand_placement_density_merge_data_plot")}
if (!dir.exists("ligand_placement_density_merge_plot")) {dir.create("ligand_placement_density_merge_plot")}
for (j in 1:nrow(df_structure_RMSD)) {
  pdb<-read.pdb(paste0("str_fin/",df_structure_RMSD$name.x[j]))  
  df_pdb<-pdb$atom
  
  df_pdb<-df_pdb%>%mutate(x=round(x,digits = 0))
  df_pdb<-df_pdb%>%mutate(y=round(y,digits = 0))
  df_pdb<-df_pdb%>%mutate(z=round(z,digits = 0))
  df_pdb<-df_pdb%>%select(x,y,z)
  df_pdb<-unique(df_pdb)
  
  write.csv(df_pdb,paste0("ligand_placement_density/",df_structure_RMSD$name.x[j],".csv"),row.names = F)
  
}
df_structure_RMSD_tost<-df_structure_RMSD%>%select(receptor,ligand)
df_structure_RMSD_tost<-unique(df_structure_RMSD_tost)
for (j in 1:nrow(df_structure_RMSD_tost)) {
  df_structure_RMSD_part<-df_structure_RMSD%>%filter(receptor==df_structure_RMSD_tost$receptor[j])%>%
    filter(ligand==df_structure_RMSD_tost$ligand[j])
  df_density<-read.csv(paste0("ligand_placement_density/",
                              df_structure_RMSD$name.x[1],".csv"),stringsAsFactors = F)
  for (p in 2:nrow(df_structure_RMSD_part)) {
    df_density_add<-read.csv(paste0("ligand_placement_density/",
                                    df_structure_RMSD$name.x[p],".csv"),stringsAsFactors = F)
    df_density<-rbind(df_density,df_density_add)
    df_density<-unique(df_density)
  }

#          df_pore_add<-df_pore[!df_pore$z%in%df_pore_size$z,]
  write.csv(df_density,paste0("ligand_placement_density_merge/",df_structure_RMSD_tost$receptor[j],"_",
                              df_structure_RMSD_tost$ligand[j],".csv"),row.names = F)
  
}
j<-1
v_sep<-seq(from=-60,to=60,by=10)

j
for (j in 1:nrow(df_structure_RMSD_tost)){
  df_density<-read.csv(paste0("ligand_placement_density_merge/",df_structure_RMSD_tost$receptor[j],"_",
                              df_structure_RMSD_tost$ligand[j],".csv"),stringsAsFactors = F)
  
  df_density<-df_density%>%group_by(z)%>%mutate(size=n())
  df_density<-df_density%>%select(z,size)
  df_pore_size<-unique(df_density)
  df_pore_add<-df_pore[!df_pore$z%in%df_pore_size$z,]
  df_pore_add<-df_pore_add%>%filter(z<max(df_pore_size$z))
  df_pore_add<-df_pore_add%>%filter(z>min(df_pore_size$z))
  df_pore_size<-rbind(df_pore_size,df_pore_add)
  df_pore_size<-df_pore_size%>%mutate(receptor=df_structure_RMSD_tost$receptor[j])
  df_pore_size<-df_pore_size%>%mutate(ligand=df_structure_RMSD_tost$ligand[j])

  write.csv(df_pore_size,paste0("ligand_placement_density_merge_data_plot/",df_structure_RMSD_tost$receptor[j],"_",
                              df_structure_RMSD_tost$ligand[j],".csv"),row.names = F)
}
j<-1
df_density<-read.csv(paste0("ligand_placement_density_merge_data_plot/",df_structure_RMSD_tost$receptor[j],"_",
                            df_structure_RMSD_tost$ligand[j],".csv"),stringsAsFactors = F)
for (j in 2:nrow(df_structure_RMSD_tost)) {
  df_density_add<-read.csv(paste0("ligand_placement_density_merge_data_plot/",df_structure_RMSD_tost$receptor[j],"_",
                                  df_structure_RMSD_tost$ligand[j],".csv"),stringsAsFactors = F)
  df_density<-rbind(df_density,df_density_add)
  
}  


#df_density<-df_density%>%select(z,average_size,receptor,ligand)
df_density<-unique(df_density)
df_density<-df_density%>%mutate(colour="bonding")
df_density$colour[df_density$size==0]<-"not bonding"
p_test<-ggplot(data=df_density)+
  geom_line(aes(x=z,y=size))+
  geom_point(aes(x=z,y=size,colour = colour))+
  geom_hline(yintercept=0)+
  facet_grid(ligand~receptor)+
  scale_x_continuous(breaks = v_sep,labels = v_sep)+
  theme_bw()
ggsave(p_test,filename = paste0("average_pore_size.png"), width = 24, height = 36, units = c("cm"), dpi = 1000 ) 


