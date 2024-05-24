part_start = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
library("FactoMineR")
library("factoextra")
library(cluster)
#install.packages( "clustertend")
#install.packages("cluster")
library(ggpmisc)
library(cowplot)
library(clValid)
library(kohonen)
setwd(part_start)
df_all_systems<-read.csv(paste0("start/all_systems.csv"),stringsAsFactors = F)
df_all_systems<-df_all_systems%>%mutate(system_name=paste0("charmm-gui-",system_name))
part_name<-paste0(part_start,"MD_analysis/")
if(!dir.exists("MD_analysis/statistic_plot")){dir.create("MD_analysis/statistic_plot")}
if(!dir.exists("MD_analysis/cluster_data/")){dir.create("MD_analysis/cluster_data/")}
if(!dir.exists("MD_analysis/cluster_data/claster_model/")){dir.create("MD_analysis/cluster_data/claster_model/")}
if(!dir.exists("MD_analysis/cluster_data/claster_pdb/")){dir.create("MD_analysis/cluster_data/claster_pdb/")}
#if(!dir.exists("MD_analysis/cluster_data/claster_data/")){dir.create("MD_analysis/cluster_data/claster_data/")}
if(!dir.exists("MD_analysis/cluster_data/cluster_data/")){dir.create("MD_analysis/cluster_data/cluster_data/")}
if(!dir.exists("MD_analysis/cluster_data/frame_statisitc/")){dir.create("MD_analysis/cluster_data/frame_statisitc/")}

if(!dir.exists("MD_analysis/cluster_plots/")){dir.create("MD_analysis/cluster_plots/")}
if(!dir.exists("MD_analysis/cluster_plots/frame_statisitc/")){dir.create("MD_analysis/fin_plots/cluster_statisitc/")}
#fin_plots

setwd(part_name)
#df_docking_interactions<-read.csv(paste0("docking_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
#nrow(df_frame_data)
library("factoextra")
j<-1
for (j in 1:nrow(df_all_systems)) {
  df_frame_data<-read.csv(paste0("fin_data/frame_data/",df_all_systems$system_name[j],".csv"),stringsAsFactors = F)
  row.names(df_frame_data)<-df_frame_data$frame_number
  
  df<-df_frame_data%>%select(protein_Total,protein_water_Total, protein_lipid_Total)
  df<-df%>%mutate(inter_protein_lipid_Total=protein_Total-protein_water_Total- protein_lipid_Total)
  df <- scale(df)
  # Compute clValid
  library(clValid)
  clmethods <- c("hierarchical")
#  clmethods <- c("hierarchical","kmeans","pam")
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
  if(nrow(df_hierarchical)>0  ){
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

    df_best_cluster_model<-df_best_cluster_model%>%mutate(quality=sqrt((protein_Total.x-protein_Total.y)^2+
                                                                         (protein_water_Total.x-protein_water_Total.y)^2+
                                                                         (protein_lipid_Total.x-protein_lipid_Total.y)^2+
                                                                         (inter_protein_lipid_Total.x-inter_protein_lipid_Total.y)^2))
    write.csv(df_best_cluster_model,paste0("cluster_data/claster_model/cluster_",df_all_systems$system_name[j],".csv"),row.names = F)
    
    df_best_cluster_model<-df_best_cluster_model%>%group_by(cluster)%>%filter(quality==min(quality))
    write.csv(df_best_cluster_model,paste0("cluster_data/frame_statisitc/cluster_",df_all_systems$system_name[j],".csv"),row.names = F)
    
    for (i in 1:nrow(df_best_cluster_model)) {
      pdb<-read.pdb(paste0(part_start,"MD/",df_all_systems$system_name[j],"/din/pdb_second/8/",df_best_cluster_model$name[i],".pdb"))
      write.pdb(pdb, paste0(part_start,"MD_analysis/cluster_data/claster_pdb/",df_all_systems$system_name[j],"_",df_best_cluster_model$cluster[i],".pdb"))
    }
    df_best_cluster_model<-read.csv(paste0("cluster_data/frame_statisitc/cluster_",df_all_systems$system_name[j],".csv"),stringsAsFactors = F)
    df_best_cluster_model<-df_best_cluster_model%>%select(cluster,name,quality)
    df_best_cluster_model<-df_best_cluster_model%>%mutate(RMSD=NA)
    df_best_cluster_compare<-left_join(df_best_cluster_model,df_best_cluster_model,by="RMSD")
    df_best_cluster_compare<-df_best_cluster_compare%>%filter(cluster.x<cluster.y)
    i<-1
    for (i in 1:nrow(df_best_cluster_compare)) {
      pdb_1<-read.pdb(paste0(part_name,"/cluster_data/claster_pdb/",
                             df_all_systems$system_name[j],"_",df_best_cluster_compare$cluster.x[i],".pdb"))
      pdb_2<-read.pdb(paste0(part_name,"cluster_data/claster_pdb/",
                             df_all_systems$system_name[j],"_",df_best_cluster_compare$cluster.y[i],".pdb"))  
      pdb_1.inds<-atom.select(pdb_1,"calpha")
      pdb_2.inds<-atom.select(pdb_2,"calpha")
      df_best_cluster_compare$RMSD[i]<-rmsd(a = pdb_1,b=pdb_2,a.inds = pdb_1.inds,b.inds = pdb_2.inds,fit = T)
    }
    write.csv(df_best_cluster_compare,paste0(part_start,"MD_analysis/cluster_data_2/claster_compare/RMSD_",df_all_systems$system_name[j],".csv"),row.names = F)
  }
}
#    part<-paste0(part_start,'MD/',df_all_systems$system_name[j],"/")
#    for (i in 1:nrow(df_best_cluster_compare)) {
#      df_tcl<-data.frame(matrix(nrow = 1,ncol = 12))
#      df_tcl[1,1]<-paste('cd', part,'\nmol new {namd/step5_input.psf} type {psf}')
#      df_tcl[1,2]<-paste0('mol addfile {din/pdb_second/hbond_8/',df_best_cluster_compare$name.x[i],'.pdb} type {pdb} waitfor all')
#      df_tcl[1,3]<-paste0('mol addfile {din/pdb_second/hbond_8/',df_best_cluster_compare$name.y[i],'.pdb} type {pdb} waitfor all')
#      df_tcl[1,4]<-paste0('set protein [atomselect top "protein and name CA"]')
#      df_tcl[1,5]<-paste0('set n [molinfo top get numframes]')
#      df_tcl[1,6]<-paste0('set output [open ',part_name,'cluster_data/cluster_data/RMSF_',df_all_systems$system_name[j],'_',
#                          df_best_cluster_compare$cluster.x[i],'_',df_best_cluster_compare$cluster.y[i],'.txt w] ')
#      df_tcl[1,7]<-paste0('set rmsf [measure rmsf $protein]')
#      df_tcl[1,8]<-paste0('foreach x $rmsf {')
#      df_tcl[1,9]<-paste0('puts $output $x')
#      df_tcl[1,10]<-paste0('}')
#      df_tcl[1,11]<-paste0('puts "output file: $n ',part_name,'cluster_data/cluster_data/RMSF_',df_all_systems$system_name[j],'_',
#                           df_best_cluster_compare$cluster.x[i],'_',df_best_cluster_compare$cluster.y[i],'.txt"')
#      df_tcl[1,12]<-paste0('close $output')
#      df_tcl[1,13]<-paste0('mol delete all\n\n\n exit now')
#      write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/compare_clusters_RMSF_',df_all_systems$system_name[j],'_',
#                                      df_best_cluster_compare$cluster.x[i],'_',df_best_cluster_compare$cluster.y[i],'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
#      system(command = paste0("vmd -dispdev text -e ",part_start,'MD_analysis/tcl/compare_clusters_RMSF_',df_all_systems$system_name[j],'_',
#                              df_best_cluster_compare$cluster.x[i],'_',df_best_cluster_compare$cluster.y[i],'.tcl'),ignore.stdout=T,wait = T)
#    }
#    df_RMSF <- read.table(paste0(part_name,'cluster_data/cluster_data/RMSF_',df_all_systems$system_name[j],'_',
#                                 df_best_cluster_compare$cluster.x[1],'_',df_best_cluster_compare$cluster.y[1],'.txt'), sep="", header=F, na.strings ="", stringsAsFactors= F)
#    colnames(df_RMSF)<-"RMSF"
#    df_RMSF <- df_RMSF%>% mutate(resid=c(1:nrow(df_RMSF)))
#    df_RMSF<-df_RMSF%>%mutate(cluster_comparition=paste(df_best_cluster_compare$cluster.x[1],'-',df_best_cluster_compare$cluster.y[1]))
#    
#    if(nrow(df_best_cluster_compare)>1){
#      for (i in 2:nrow(df_best_cluster_compare)) {    
#        df_RMSF_add <- read.table(paste0(part_name,'cluster_data/cluster_data/RMSF_',df_all_systems$system_name[j],'_',
#                                         df_best_cluster_compare$cluster.x[i],'_',df_best_cluster_compare$cluster.y[i],'.txt'), sep="", header=F, na.strings ="", stringsAsFactors= F)
#        colnames(df_RMSF_add)<-"RMSF"
#        df_RMSF_add <- df_RMSF_add%>% mutate(resid=c(1:nrow(df_RMSF)))
#        df_RMSF_add<-df_RMSF_add%>%mutate(cluster_comparition=paste(df_best_cluster_compare$cluster.x[i],'-',df_best_cluster_compare$cluster.y[i]))
#        df_RMSF<-rbind(df_RMSF,df_RMSF_add)
#      }
#    }
#    df_RMSD<-df_best_cluster_compare%>%select(cluster.x,cluster.y,RMSD)
#    colnames(df_RMSD)<-c("cluster 1","cluster 2","RMSD, A")
#    p_RMSD<-ggplot()+ labs(x="",y="")+
#      annotate(geom = "table",
#               x = 9,
#               y = 3,
#               label = list(df_RMSD))+theme_minimal()+
#      scale_x_continuous(breaks = NULL,labels = NULL)+
#      scale_y_continuous(breaks = NULL,labels = NULL)
#    df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors =  F)
#    df_topology<-df_topology%>%filter(topology=='TMD')
#    v_sep<-max(round(df_RMSF$resid,digits = -1))
#    v_sep<-seq(from=0,to=v_sep,by=20)
#    if(nrow(df_best_cluster_compare)>1){
#      p_RMSF<-ggplot()+   
#        geom_rect(aes(xmin=seq_beg-0.5,xmax=seq_end+0.5,ymin=-Inf,ymax=Inf),alpha=0.5,data=df_topology)+
#        geom_line(data=df_RMSF,aes(x=resid,y=RMSF,colour=cluster_comparition))+
#        scale_x_continuous(breaks = v_sep,labels = v_sep)+
#        theme_bw()+theme(legend.position = "bottom")
#    }else{
#      p_RMSF<-ggplot()+   
#        geom_rect(aes(xmin=seq_beg-0.5,xmax=seq_end+0.5,ymin=-Inf,ymax=Inf),alpha=0.5,data=df_topology)+
#        geom_line(data=df_RMSF,aes(x=resid,y=RMSF))+
#        scale_x_continuous(breaks = v_sep,labels = v_sep)+
#        theme_bw()+
#        scale_y_continuous(breaks=c(0:10),labels=c(0:10),limits=c(0,10))
#    }
#    p_RMSF<-ggdraw() +
#      draw_plot(p_RMSF)+
#      draw_plot(p_RMSD,
#                x = 0.55, y = 0.35, hjust = 0, vjust = 0, 
#                halign = 0.8, valign = 0, width = 0.4)
#    
#    # Visualize the tree
#    p_tree<-fviz_dend(res.hk, cex = 0.6, palette = "jco",
##                      rect = TRUE, rect_border = "jco", rect_fill = TRUE,
#                      geom = c("point"), show_labels = F)
#    p_annotate<- ggplot()+ labs(x="",y="")+
#      annotate(geom = "table",
#               x = 9,
#               y = 3,
#               label = list(df_cluster_statistic))+theme_minimal()+
#      scale_x_continuous(breaks = NULL,labels = NULL)+
#      scale_y_continuous(breaks = NULL,labels = NULL)
#    p_tree_merge<-ggdraw() +
#      draw_plot(p_tree)+
#      draw_plot(p_annotate,
#                x = 0.6, y = 0.35, hjust = 0, vjust = 0, 
#                halign = 0.8, valign = 0, width = 0.4)
    #    p_merge
    #                    ggtheme = theme_cleveland())
    # Visualize the hkmeans final clusters
#    p_clust<-fviz_cluster(res.hk, palette = "jco", repel = TRUE,
#                          geom = c("point"),
#                          ggtheme = theme_bw())
#    p_merge<-plot_grid(p_clust,p_tree_merge,p_RMSF,rel_widths=c(1,1,2),nrow=1)
#    ggsave(p_merge,file=paste0("cluster_plots/frame_statisitc/cluster_quantity_kmeans_",df_all_systems$Structure[j],"_",
#                               df_all_systems$Membrane[j],".png"), width = 60, height = 10, units = c("cm"), dpi = 200 )
#  }else {print(j)}
#}
