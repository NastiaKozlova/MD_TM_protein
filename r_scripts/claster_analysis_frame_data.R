part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
library("FactoMineR")
library("factoextra")
library(ggpmisc)
library(cowplot)
setwd(part_start)
df_all_systems<-read.csv(paste0("start/all_systems.csv"),stringsAsFactors = F)
df_all_systems<-df_all_systems%>%mutate(system_name=paste0("charmm-gui-",system_name))
part_name<-paste0(part_start,"MD_analysis/fin_data/")
if(!dir.exists("MD_analysis/statistic_plot")){dir.create("MD_analysis/statistic_plot")}
if(!dir.exists("MD_analysis/fin_data/claster_model/")){dir.create("MD_analysis/fin_data/claster_model/")}
if(!dir.exists("MD_analysis/fin_data/claster_pdb/")){dir.create("MD_analysis/fin_data/claster_pdb/")}

setwd(part_name)
#df_docking_interactions<-read.csv(paste0("docking_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
#nrow(df_frame_data)
j<-2

for (j in 1:nrow(df_all_systems)) {
  #prepare data for analysis
  df_frame_data<-read.csv(paste0("frame_data/",df_all_systems$system_name[j],".csv"),stringsAsFactors = F)
  row.names(df_frame_data)<-df_frame_data$number
  df_frame_data<-df_frame_data%>%select(ramachadran, protein_Total, protein_water_Total,
                                        protein_lipid_Total,SASA_protein,RMSD_protein)
  res.pca <- PCA(df_frame_data, scale.unit = TRUE, ncp = 5, graph = F)
  res.hcpc <- HCPC(res.pca, graph = FALSE, nb.clust = 0, min = 3, max = NULL)
  #analysis of clasters
  data<-res.hcpc[["data.clust"]]
  data<-data%>%mutate(frame=row.names(data))
  df_data<-as.data.frame(data)
  df_data<-df_data%>%group_by(clust)%>%mutate(size=n())
  df_data<-ungroup(df_data)
  #select the best model in the claster 
  df_best_cluster_model<-df_data%>%group_by(clust)%>%mutate(min_ramachandran=min(ramachadran))
  df_best_cluster_model<-df_best_cluster_model%>%filter(min_ramachandran==ramachadran)
  
  df_best_cluster_model<-df_best_cluster_model%>%group_by(clust)%>%mutate(min_protein_Total=min(protein_Total))
  df_best_cluster_model<-df_best_cluster_model%>%filter(min_protein_Total==protein_Total)
  df_best_cluster_model<-ungroup(df_best_cluster_model)
  write.csv(df_best_cluster_model,paste0(part_start,"MD_analysis/fin_data/claster_model/",df_all_systems$system_name[j],".csv"),row.names=F)

  # describe PCA
  eig.val <- get_eigenvalue(res.pca)
  #PCA analysis and making plot
  p_pca_pa<-fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))+theme_bw()
  
  # Contributions to the dimensions
  var <- get_pca_var(res.pca)
  df_contrib<-head(var$contrib)
  a<-row.names(df_contrib)
  df_contrib<-as.data.frame(df_contrib)
  df_contrib_add<-data.frame(matrix(nrow = nrow(df_contrib),ncol = 1))
  colnames(df_contrib_add)<-c("parameters")
  df_contrib_add$parameters<-a
  df_contrib<-cbind(df_contrib_add,df_contrib)
  #row.names(df_contrib)<-a
  p_contrib<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_contrib))+theme_bw()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  df_best_cluster<-df_best_cluster_model%>%select(clust,size)
  p_stat<-plot_grid(p_pca_pa,p_contrib,ncol=2)

  p_annotate<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_best_cluster))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p<-fviz_dend(res.hcpc, show_labels = F)
  p_clustering<- ggdraw(p) + 
    draw_plot(p_annotate, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
  p_fin<-plot_grid(p_stat,p_clustering,ncol = 1,rel_heights=c(0.5,1))
  ggsave(p_fin,filename = paste0(part_start,"MD_analysis/statistic_plot/",df_all_systems$system_name[j],".png"),width=15,height=15)
  
}
i<-1
for (j in 1:nrow(df_all_systems)) {
  df_best_cluster_model<-read.csv(paste0(part_start,"MD_analysis/fin_data/claster_model/",df_all_systems$system_name[j],".csv"),
                                  stringsAsFactors = F)
  for (i in 1:nrow(df_best_cluster_model)) {
    pdb<-read.pdb(paste0(part_start,"MD/",df_all_systems$system_name[j],"/din/pdb_second/8/frame_",df_best_cluster_model$frame[i],".pdb"))
    write.pdb(pdb, paste0(part_start,"MD_analysis/fin_data/claster_pdb/",df_all_systems$system_name[j],"_",df_best_cluster_model$clust[i],".pdb"))
  }
}
for (j in 1:nrow(df_all_systems)) {
  df_best_cluster_model<-read.csv(paste0(part_start,"MD_analysis/fin_data/claster_model/",df_all_systems$system_name[j],".csv"),
                                  stringsAsFactors = F)
  df_best_cluster_model<-df_best_cluster_model%>%select(clust,frame)
  df_best_cluster_model<-df_best_cluster_model%>%mutate(RMSD=NA)
  
  df_best_cluster<-left_join(df_best_cluster_model,df_best_cluster_model,by="RMSD")
  for (i in 1:nrow(df_best_cluster)) {
    pdb1<-read.pdb(paste0(part_start,"MD/",df_all_systems$system_name[j],"/din/pdb_second/8/frame_",df_best_cluster$frame.x[i],".pdb"))
    pdb2<-read.pdb(paste0(part_start,"MD/",df_all_systems$system_name[j],"/din/pdb_second/8/frame_",df_best_cluster$frame.y[i],".pdb"))
    df_best_cluster$RMSD[i]<-rmsd(pdb1,pdb2)
  }
}
  