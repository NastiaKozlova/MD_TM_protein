part_start = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
library("FactoMineR")
library("factoextra")
library(cluster)
#install.packages( "clustertend")
#install.packages("cluster")
library(mclust)
library(ggpmisc)
library(cowplot)
library(clValid)
library(kohonen)
setwd(part_start)
df_all_systems<-read.csv(paste0("start/all_systems.csv"),stringsAsFactors = F)
df_all_systems<-df_all_systems%>%mutate(system_name=paste0("charmm-gui-",system_name))
part_name<-paste0(part_start,"MD_analysis/")
#if(!dir.exists("MD_analysis/statistic_plot")){dir.create("MD_analysis/statistic_plot")}
if(!dir.exists("MD_analysis/cluster_data_all/")){dir.create("MD_analysis/cluster_data_all/")}
if(!dir.exists("MD_analysis/cluster_data_all/claster_model/")){dir.create("MD_analysis/cluster_data_all/claster_model/")}
if(!dir.exists("MD_analysis/cluster_data_all/claster_pdb/")){dir.create("MD_analysis/cluster_data_all/claster_pdb/")}
if(!dir.exists("MD_analysis/cluster_data_all/claster_compare/")){dir.create("MD_analysis/cluster_data_all/claster_compare/")}
if(!dir.exists("MD_analysis/cluster_data_all/cluster_data/")){dir.create("MD_analysis/cluster_data_all/cluster_data/")}
if(!dir.exists("MD_analysis/cluster_data_all/frame_statisitc/")){dir.create("MD_analysis/cluster_data_all/frame_statisitc/")}

if(!dir.exists("MD_analysis/cluster_plots_all/")){dir.create("MD_analysis/cluster_plots_all/")}
if(!dir.exists("MD_analysis/cluster_plots_all/cluster_statisitc/")){dir.create("MD_analysis/cluster_plots_all/cluster_statisitc/")}
#fin_plots

setwd(part_name)
#df_docking_interactions<-read.csv(paste0("docking_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
#nrow(df_frame_data)
library("factoextra")
j<-1
df_frame_data<-read.csv(paste0("fin_data/frame_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
df_frame_data<-df_frame_data%>%mutate(system_name=df_all_systems$system_name[1])
for (j in 2:nrow(df_all_systems)) {
  df_frame_data_add<-read.csv(paste0("fin_data/frame_data/",df_all_systems$system_name[j],".csv"),stringsAsFactors = F)
  df_frame_data_add<-df_frame_data_add%>%mutate(system_name=df_all_systems$system_name[j])
  df_frame_data<-rbind(df_frame_data,df_frame_data_add)
}
df_frame_data<-df_frame_data%>%mutate(name=paste0(system_name,"_",frame_number))
row.names(df_frame_data)<-df_frame_data$name
  
df<-df_frame_data%>%select( ramachadran,  protein_Total,  
                            protein_water_Total, protein_lipid_Total,
                            SASA_protein, RMSD_protein)
df<-df%>%mutate(inter_protein_Total=protein_Total-protein_water_Total- protein_lipid_Total)
df <- scale(df)
  
res.pca<-PCA(df,graph = F)
eig.val<-get_eigenvalue(res.pca)
fviz_eig(res.pca,addlabels = T,ylim=c(0,50))
fviz_pca_var(res.pca,col.var = "black")
var<-get_pca_var(res.pca = res.pca)
  
library("corrplot")
corrplot(var$cos2)



  # Compute clValid
  library(clValid)
  #clmethods <- c("hierarchical")
  clmethods <- c("hierarchical","kmeans","pam")
  intern_1 <- clValid(df, nClust = 2:10,
                      clMethods = clmethods, 
                      metric = "euclidean",
                      validation =c("internal","stability"),
                      maxitems=nrow(df),verbose=F)
  # Summary
  a<- optimalScores(intern_1)
  df_a<-as.data.frame(a)
  df_summary<-df_a%>%group_by(Method,Clusters)%>%summarise(quantity=(n()/nrow(df_a)))
  df_hierarchical<-df_summary%>%filter(Method=="hierarchical")
#  if(nrow(df_hierarchical)>0  ){
    df_hierarchical<-df_hierarchical%>%filter(quantity==max(quantity))
    # Compute hierarchical k-means clustering
    library(factoextra)
    res.hk <-hkmeans(df, df_hierarchical$Clusters[1])
    # Elements returned by hkmeans()
    #names(res.hk)
    
    df_all<-as.data.frame(res.hk$cluster)
    colnames(df_all)<-"cluster"
    
    df_all<-df_all%>%group_by(cluster)%>%mutate(size=round(n()/nrow(df_all)*100,digits = 2))%>%arrange(cluster)
    df_cluster_statistic<-unique(df_all)
    df_center<-as.data.frame(res.hk$centers)
    df_center<-df_center%>%mutate(cluster=rownames(df_center))
    
    df_all<-as.data.frame(res.hk$cluster)
    colnames(df_all)<-"cluster"
    df_all <-df_all%>%mutate(name=row.names(df_all))
    df_temp<-as.data.frame(df)
    df_temp<-df_temp %>%mutate(name=row.names(df_temp))
    df_all<-left_join(df_all,df_temp,by="name")
    df_all<-df_all%>%mutate(cluster=as.character(cluster))
    df_best_cluster_model<-left_join(df_all,df_center,by="cluster")
    
    df_best_cluster_model<-df_best_cluster_model%>%mutate(quality=sqrt((ramachadran.x-ramachadran.y)^2+
                                                                         (protein_Total.x-protein_Total.y)^2+
                                                                         (protein_water_Total.x-protein_water_Total.y)^2+
                                                                         (protein_lipid_Total.x-protein_lipid_Total.y)^2+
                                                                         (inter_protein_Total.x-inter_protein_Total.y)^2))
    write.csv(df_best_cluster_model,paste0("cluster_data_all/claster_model/cluster_",df_all_systems$system_name[j],".csv"),row.names = F)
    
    df_best_cluster_model<-df_best_cluster_model%>%group_by(cluster)%>%filter(quality==min(quality))
    write.csv(df_best_cluster_model,paste0("cluster_data_all/frame_statisitc/cluster_",df_all_systems$system_name[j],".csv"),row.names = F)
    
    for (i in 1:nrow(df_best_cluster_model)) {
      pdb<-read.pdb(paste0(part_start,"MD/",df_all_systems$system_name[j],"/din/pdb_second/8/",df_best_cluster_model$name[i],".pdb"))
      write.pdb(pdb, paste0(part_start,"MD_analysis/cluster_data_all/claster_pdb/",df_all_systems$system_name[j],"_",df_best_cluster_model$cluster[i],".pdb"))
    }
    df_best_cluster_model<-read.csv(paste0("cluster_data_all/frame_statisitc/cluster_",df_all_systems$system_name[j],".csv"),stringsAsFactors = F)
    df_best_cluster_model<-df_best_cluster_model%>%select(cluster,name,quality)
    df_best_cluster_model<-df_best_cluster_model%>%mutate(RMSD=NA)
    df_best_cluster_compare<-left_join(df_best_cluster_model,df_best_cluster_model,by="RMSD")
    df_best_cluster_compare<-df_best_cluster_compare%>%filter(cluster.x<cluster.y)
    i<-1
    for (i in 1:nrow(df_best_cluster_compare)) {
      pdb_1<-read.pdb(paste0(part_name,"/cluster_data_all/claster_pdb/",
                             df_all_systems$system_name[j],"_",df_best_cluster_compare$cluster.x[i],".pdb"))
      pdb_2<-read.pdb(paste0(part_name,"cluster_data_all/claster_pdb/",
                             df_all_systems$system_name[j],"_",df_best_cluster_compare$cluster.y[i],".pdb"))  
      pdb_1.inds<-atom.select(pdb_1,"calpha")
      pdb_2.inds<-atom.select(pdb_2,"calpha")
      df_best_cluster_compare$RMSD[i]<-rmsd(a = pdb_1,b=pdb_2,a.inds = pdb_1.inds,b.inds = pdb_2.inds,fit = T)
      #    v_rmsf<-rmsf(pdb_1)
      pdb<-read.pdb(paste0(part_start,"MD/",df_all_systems$system_name[j],"/din/pdb_second/8/",df_best_cluster_model$name[i],".pdb"))
      write.pdb(pdb, paste0(part_start,"MD_analysis/cluster_data_all/claster_pdb/",df_all_systems$system_name[j],"_",df_best_cluster_model$cluster[i],".pdb"))
    }
    write.csv(df_best_cluster_compare,paste0(part_start,"MD_analysis/cluster_data_all/claster_compare/RMSD_",df_all_systems$system_name[j],".csv"),row.names = F)
  }
}
