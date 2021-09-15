part_start <- commandArgs(trailingOnly=TRUE)
library(dplyr)
library(ggplot2)
library(cowplot)
library(bio3d)
library(bio3d)
library(ggplot2)
library(dplyr)
test_10<-seq(from=0,to=1000,by=10)

part<-paste0(part_start,"MD_analysis/")
parta<-paste0(part_start,"MD/")
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)

df_all_systems$system_name<-as.character(df_all_systems$system_name)
#df_all_systems<-df_all_systems%>%mutate(Progress=="notdone")
#df_all_systems_done<-df_all_systems%>%filter(Progress=="done")
#df_all_systems<-rbind(df_all_systems_productive_run,df_all_systems_done)
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
main_part<-c(8)
i<-nrow(df_all_systems)
main<-main_part[2]
if (!dir.exists(paste0(part,"fin_data"))) {dir.create(paste0(part,"fin_data"))}
if (!dir.exists(paste0(part,"fin_data/frame_data"))) {dir.create(paste0(part,"fin_data/frame_data"))}
if (!dir.exists(paste0(part,"fin_data/str_data"))) {dir.create(paste0(part,"fin_data/str_data"))}
if (!dir.exists(paste0(part,"fin_data/str_data_sep"))) {dir.create(paste0(part,"fin_data/str_data_sep"))}
if (!dir.exists(paste0(part,"fin_plots"))) {dir.create(paste0(part,"fin_plots"))}
if (!dir.exists(paste0(part,"fin_plots/test_docking_plots/"))) {dir.create(paste0(part,"fin_plots/test_docking_plots/"))}
if (!dir.exists(paste0(part,"fin_plots/str_plots"))) {dir.create(paste0(part,"fin_plots/str_plots"))}
if (!dir.exists(paste0(part,"fin_plots/str_XYZ"))) {dir.create(paste0(part,"fin_plots/str_XYZ"))}
if (!dir.exists(paste0(part,"fin_plots/frame_plots_histo"))) {dir.create(paste0(part,"fin_plots/frame_plots_histo"))}
if (!dir.exists(paste0(part,"fin_plots/docking_plots"))) {dir.create(paste0(part,"fin_plots/docking_plots"))}
#if (!dir.exists(paste0(part,"fin_plots/lenght_plots"))) {dir.create(paste0(part,"fin_plots/lenght_plots"))}
i<-1
main<-main_part[1]
for (main in main_part) {
  df_fin<-read.csv(paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[1],"_",main,".csv"),stringsAsFactors = F)
  df_fin<-df_fin%>%mutate(main=main)
  df_fin<-df_fin%>%mutate(system=df_all_systems$fin_name[1])
  for (i in 2:nrow(df_all_systems)) {
    df_fin_add<-read.csv(paste0(part,"fin_data/frame_data/",df_all_systems$fin_name[i],"_",main,".csv"),stringsAsFactors = F)
    df_fin_add<-df_fin_add%>%mutate(main=main)
    df_fin_add<-df_fin_add%>%mutate(system=df_all_systems$fin_name[i])
    df_fin<-rbind(df_fin,df_fin_add)
  }
}
df_fin<-df_fin%>%mutate(protein_protein_Total=protein_Total-protein_water_Total-protein_lipid_Total)
df_fin<-left_join(df_fin,df_all_systems,by=c("system"="fin_name"))
#df_fin<-df_fin%>%filter(frame>10)
df_fin<-df_fin%>%mutate(system=paste0(Structure,"_",Membrane))
p_protein_lipid_histo<-ggplot(data = df_fin)+
  labs(title=paste("protein_lipid"), x = "protein_lipid (kcal/mol)") +theme_bw()+
  geom_freqpoly(aes(x =protein_lipid_Total,colour=system))+guides(colour = "none")
p_protein_water_histo<-ggplot(data = df_fin)+
  labs(title=paste("protein_water"), x = "protein_water (kcal/mol)") +theme_bw()+
  geom_freqpoly(aes(x =protein_water_Total,colour=system))+guides(colour = "none")
p_protein_protein_histo<-ggplot(data = df_fin)+
  labs(title=paste("protein_protein"), x = "protein_protein (kcal/mol)") +theme_bw()+
  geom_freqpoly(aes(x =protein_protein_Total,colour=system))+guides(colour = "none")
p_protein_histo<-ggplot(data = df_fin)+
  labs(title=paste("protein"), x = "protein (kcal/mol)") +theme_bw()+
  geom_freqpoly(aes(x =protein_Total,colour=system))+theme(legend.position = "bottom")
p_all<-plot_grid(p_protein_lipid_histo,p_protein_water_histo,       p_protein_protein_histo,      p_protein_histo,
                 rel_heights=c(1,1,1,1),
                 nrow=4,  labels = c("A","B","C","D"),align="hv",ncol = 1,axis="bt")


ggsave(p_all,   filename = paste0(part,"fin_plots/tost_",main,".png"), width = 20, height = 25, units = c("cm"), dpi = 200 ) 
M <- aov(protein_Total ~ system, data = df_fin)
df_protein_Total<-TukeyHSD(M)[[1]]

M <- aov(protein_protein_Total ~ system, data = df_fin)
df_protein_protein_Total<-TukeyHSD(M)[[1]]

M <- aov(protein_lipid_Total ~ system, data = df_fin)
df_protein_lipid_Total<-TukeyHSD(M)[[1]]

M <- aov(protein_water_Total ~ system, data = df_fin)
df_protein_water_Total<-TukeyHSD(M)[[1]]
