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
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
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
  labs(x = "Energy (kcal/mol)",title="Energy of protein lipid interactions",) +theme_bw()+
  geom_freqpoly(aes(x =protein_lipid_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
p_protein_water_histo<-ggplot(data = df_fin)+
  labs(x = "Energy (kcal/mol)",title="Energy of protein water interactions") +theme_bw()+
  geom_freqpoly(aes(x =protein_water_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
p_protein_protein_histo<-ggplot(data = df_fin)+
  labs(x = "Energy (kcal/mol)",title="Energy of interprotein interactions") +theme_bw()+
  geom_freqpoly(aes(x =protein_protein_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
p_protein_histo<-ggplot(data = df_fin)+
  labs(x = "Energy (kcal/mol)",title="Total protein energy") +theme_bw()+
  geom_freqpoly(aes(x =protein_Total,colour=Structure,linetype=Membrane))+theme(legend.position = "bottom")
legend_test<-get_legend(p_protein_histo)

p_all<-plot_grid(p_protein_lipid_histo+theme(legend.position = "none"),
                 p_protein_water_histo+theme(legend.position = "none"),  
                 p_protein_protein_histo+theme(legend.position = "none"),
                 p_protein_histo+theme(legend.position = "none"),
                 rel_heights=c(1,1),
                 nrow=2,  labels = c("A","B","C","D"),align="hv",ncol = 2,axis="bt")
p_test<-plot_grid(p_all,legend_test,nrow=2,rel_heights=c(1,0.1))
