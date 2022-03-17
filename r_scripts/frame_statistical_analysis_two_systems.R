part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggpmisc)
library(cowplot)
library(rlang)
setwd(part_start)
df_all_systems<-read.csv(paste0("start/all_systems.csv"),stringsAsFactors = F)
df_all_systems<-df_all_systems%>%mutate(system_name=paste0("charmm-gui-",system_name))
part_name<-paste0(part_start,"MD_analysis/fin_data/")
if(!dir.exists("MD_analysis/statistic_plot")){dir.create("MD_analysis/statistic_plot")}
setwd(part_name)
df_docking_interactions<-read.csv(paste0("docking_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
df_frame_data<-read.csv(paste0("frame_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
df_frame_data<-df_frame_data%>%mutate(system_name=df_all_systems$system_name[1])
i<-2
df_aminoacids<-read.csv(paste0("str_data/",df_all_systems$system_name[1],".csv"),stringsAsFactors = F)
df_aminoacids<-df_aminoacids%>%mutate(system_name=df_all_systems$system_name[1])
#nrow(df_frame_data)
if(nrow(df_all_systems)>1){
  for (i in 2:nrow(df_all_systems)) {
    df_frame_data_add<-read.csv(paste0("frame_data/",df_all_systems$system_name[i],".csv"),stringsAsFactors = F)
    df_frame_data_add<-df_frame_data_add%>%mutate(system_name=df_all_systems$system_name[i])
    df_frame_data<-rbind(df_frame_data,df_frame_data_add)
    
    df_aminoacids_add<-read.csv(paste0("str_data/",df_all_systems$system_name[i],".csv"),stringsAsFactors = F)
    df_aminoacids_add<-df_aminoacids_add%>%mutate(system_name=df_all_systems$system_name[i])
    df_aminoacids<-rbind(df_aminoacids,df_aminoacids_add)
  }
}

df_frame_data<-left_join(df_frame_data,df_all_systems,by = "system_name")
df_frame_data<-df_frame_data%>%mutate(system=paste0(Structure,"_",Membrane))

df_aminoacids<-left_join(df_aminoacids,df_all_systems,by = "system_name")
df_aminoacids<-df_aminoacids%>%mutate(system=paste0(Structure,"_",Membrane))

if(length(unique(df_frame_data$system))==2){
  wilcox_test_protein_Total <- df_frame_data %>% 
    wilcox_test(protein_Total ~ system)%>%
    add_significance()
  df_protein_Total<-wilcox_test_protein_Total%>%select(group1,group2,p)
  p<-ggplot(data=df_frame_data)+    
    labs(title=paste0("Wilcoxon test, n = ", nrow(df_frame_data) ),
      x="Рrotein energy, kkal/mol")+
    geom_density(aes(x=protein_Total,colour=system))+
    theme_bw()+
    theme(legend.position="none")
  colnames(df_protein_Total)<-c("system 1","system 2","p-value")
  p1<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_protein_Total))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p_protein_Total<- ggdraw(p) + 
    draw_plot(p1, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
  wilcox_test_protein_water_Total <- df_frame_data %>% 
    wilcox_test(protein_water_Total ~ system)%>%
    add_significance()
  df_protein_water_Total<-wilcox_test_protein_water_Total%>%select(group1,group2,p)
  p<-ggplot(data=df_frame_data)+    
    labs(x="Energy of protein-water interactions, kkal/mol")+
    geom_density(aes(x=protein_water_Total,colour=system))+
    theme_bw()+
    theme(legend.position="none")#+
  colnames(df_protein_water_Total)<-c("system 1","system 2","p-value")
  p1<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_protein_water_Total))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p_protein_water_Total<- ggdraw(p) + 
    draw_plot(p1, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
  wilcox_test_protein_lipid_Total <- df_frame_data %>% 
    wilcox_test(protein_lipid_Total ~ system)%>%
    add_significance()
  df_protein_lipid_Total<-wilcox_test_protein_lipid_Total%>%select(group1,group2,p)
  p<-ggplot(data=df_frame_data)+    
    labs(x="Energy of protein-lipid interactions, kkal/mol")+
    geom_density(aes(x=protein_lipid_Total,colour=system))+
    theme_bw()+
    theme(legend.position="none")#+
  colnames(df_protein_lipid_Total)<-c("system 1","system 2","p-value")
  p1<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_protein_lipid_Total))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p_protein_lipid_Total<- ggdraw(p) + 
    draw_plot(p1, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
  wilcox_test_SASA_protein <- df_frame_data %>% 
    wilcox_test(SASA_protein ~ system)%>%
    add_significance()
  df_SASA_protein<-wilcox_test_SASA_protein%>%select(group1,group2,p)
  p<-ggplot(data=df_frame_data)+    
    labs(x="SASA, AA^2")+
    geom_density(aes(x=SASA_protein,colour=system))+
    theme_bw()+
    theme(legend.position="none")
  colnames(df_SASA_protein)<-c("system 1","system 2","p-value")
  p1<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_SASA_protein))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p_SASA_protein<- ggdraw(p) + 
    draw_plot(p1, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
  wilcox_test_RMSD_protein <- df_frame_data %>% 
    wilcox_test(RMSD_protein ~ system)%>%
    add_significance()
  df_RMSD_protein<-wilcox_test_RMSD_protein%>%select(group1,group2,p)
  p<-ggplot(data=df_frame_data)+    
    labs(#subtitle = get_test_label(wilcox_test_statiotion, detailed = TRUE),
      x="RMSD, AA")+
    geom_density(aes(x=RMSD_protein,colour=system))+
    theme_bw()+
    theme(legend.position="none")#+
  colnames(df_RMSD_protein)<-c("system 1","system 2","p-value")
  p1<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_RMSD_protein))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p_RMSD_protein<- ggdraw(p) + 
    draw_plot(p1, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)
  wilcox_test_ramachadran <- df_frame_data %>% 
    wilcox_test(ramachadran ~ system)%>%
    add_significance()
  df_ramachadran<-wilcox_test_ramachadran%>%select(group1,group2,p)
  p<-ggplot(data=df_frame_data)+    
    labs(#subtitle = get_test_label(wilcox_test_statiotion, detailed = TRUE),
      x="Ramachandran")+
    geom_freqpoly(aes(x=ramachadran,colour=system,after_stat(density)),binwidth=1)+
    theme_bw()+
    theme(legend.position="bottom")#+
  legend <- get_legend(
    # create some space to the left of the legend
    p +  
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  p<-p + theme(legend.position="none")
  colnames(df_ramachadran)<-c("system 1","system 2","p-value")
  p1<- ggplot()+ labs(x="",y="")+
    annotate(geom = "table",
             x = 9,
             y = 3,
             label = list(df_ramachadran))+theme_minimal()+
    scale_x_continuous(breaks = NULL,labels = NULL)+
    scale_y_continuous(breaks = NULL,labels = NULL)
  p_ramachadran<- ggdraw(p) + 
    draw_plot(p1, x = 0.75, y = 0.95, hjust = 1, vjust = 1, width = 0.13, height = 0.2)

  p_statistic<-plot_grid(p_protein_Total,
                         p_protein_water_Total,
                         p_protein_lipid_Total,
                         p_SASA_protein,
                         p_RMSD_protein,
                         p_ramachadran,
                         ncol =2)
  p_statistic<-p_statistic + draw_grob(legend,  hjust = -100, vjust = -100)
  ggsave(p_statistic,filename = paste0(part_start,"MD_analysis/statistic_plot/frame_parameres_density.png"),width =30,height = 30,dpi = 300,units = "cm")
}
