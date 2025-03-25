part_name <- commandArgs(trailingOnly=TRUE)
library(bio3d)
library(dplyr)
library(ggplot2)
#v_rmsd<-4

setwd(part_name)
setwd("din")


df_all<-read.csv(paste0("df_merge_structure_log_center.csv"),stringsAsFactors = F)
df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand,"_",center.x))

#if (dir.exists(paste0("interaction/"))) { system(command = paste0("rm -r ",part_name,"din/interaction/"))}
#if (dir.exists(paste0("interaction_TEMP/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_TEMP/"))}
if (dir.exists(paste0("interaction_center/"))) {system(command = paste0("rm -r ",part_name,"din/interaction_center/"))}

if (!dir.exists(paste0("interaction_center/"))) { dir.create(paste0("interaction_center/"))}
i<-1
j<-1
p<-1
df_structure<-df_all%>%select(name.x, receptor, ligand, surf)
df_structure<-unique(df_structure)
for (j in 1:nrow(df_structure)) {
  df_complex<-df_all%>%filter(name.x==df_structure$name.x[j])
  pdb<-read.pdb(paste0(part_name,"receptor_start/",df_all$receptor[j],".pdb"))
  
  df_pdb<-pdb$atom
  df_pdb<-df_pdb%>%filter(elety=="CA")
  df_pdb<-df_pdb%>%mutate(number_interactions=0)
  df_pdb<-df_pdb%>%mutate(tested_structure=0)
  df_pdb<-df_pdb%>%mutate(total_structure=nrow(df_complex))
  test<-nrow(df_pdb)
  w<-0
  for (p in 1:nrow(df_complex)) {
    a<-strsplit(df_complex$models.y[p],split = ".",fixed = T)[[1]][1]
    b<-strsplit(a,split = "_",fixed = T)[[1]][2]
    if(file.exists(paste0("interaction/",df_complex$receptor_ligand[p],"/",b,".csv"))){
      df_protein<-read.csv(paste0("interaction/",df_complex$receptor_ligand[p],"/",b,".csv"),
                           stringsAsFactors = F) 
      df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]<-df_pdb$number_interactions[df_pdb$resno%in%df_protein$resid]+1
      df_pdb$tested_structure<-df_pdb$tested_structure+1
      
    }else{
      w<-w+1
      print(paste0("interaction/",df_complex$receptor_ligand[p],"/",b,".csv"))
    }
  }
  df_pdb<-df_pdb%>%mutate(persent_interactions=number_interactions/tested_structure*100)
  write.csv(df_pdb,
            paste0("interaction_center/",df_structure$name.x[j],".csv"),row.names = F)
}

df_all<-df_all%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))
df_all<-df_all%>%select(name.x,receptor,ligand,receptor_ligand)
df_all<-unique(df_all)
df_all<-df_all%>%group_by(receptor_ligand)%>%mutate(number=1:n())
df_all<-df_all%>%mutate(number=as.character(number))
df_all<-ungroup(df_all)
df_pdb<-read.csv(paste0("interaction_center/",df_all$name.x[1],".csv"),stringsAsFactors = F)
df_pdb<-df_pdb%>%mutate(name.x=df_all$name.x[1])
df_pdb<-df_pdb%>%filter(persent_interactions==100)
for (j in 2:nrow(df_all)) {
  df_pdb_add<-read.csv(paste0("interaction_center/",df_all$name.x[j],".csv"),stringsAsFactors = F)
  df_pdb_add<-df_pdb_add%>%mutate(name.x=df_all$name.x[j])
  df_pdb<-rbind(df_pdb,df_pdb_add)
  df_pdb<-df_pdb%>%filter(persent_interactions==100)
}

df_pdb<-left_join(df_pdb,df_all,by="name.x")
part_start<-strsplit(part_name,split = "/",fixed = T)[[1]]
part_start<-paste0(part_start[1:(length(part_start)-3)],collapse = "/")
part_start<-paste0(part_start,"/")
df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors = F)
v_seq<-seq(from=0,to =max(df_topology$seq_end),by=10)
p<-ggplot()+
  geom_rect(aes(xmin = seq_beg-0.5, xmax = seq_end+0.5, ymin = -Inf, ymax = Inf,fill=topology,alpha=0.1),data=df_topology)+
  geom_point(aes(x=resno,y=number),data=df_pdb)+
  theme_bw()+facet_grid(receptor_ligand~., scales = "free")+
  scale_x_continuous(breaks = v_seq,labels = v_seq)+
  guides(alpha = "none")
ggsave(p,   filename = paste0("center_interactions.png"), width = 60, height = 30, units = c("cm"), dpi = 200 ) 
