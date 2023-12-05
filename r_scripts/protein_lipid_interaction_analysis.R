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
part_name<-paste0(part_start,"MD_analysis/din/")
v_systems<-list.files(paste0(part_start,"MD/"))
#if(!dir.exists("MD_analysis/din/protein_lipid_interactions")){dir.create("MD_analysis/statistic_plot")}
setwd(part_name)

#library("factoextra")
j<-1
for (j in 1:length(v_systems)) {
    
    df_frame_data<-read.csv(paste0(part_name,v_systems[j],"/Energy_protein_lipid.csv"),stringsAsFactors = F)
    
    #df_frame_data<-df_frame_data%>%filter(Elec<0)
    #df_frame_data<-df_frame_data%>%filter(VdW<0)
    #df_frame_data<-df_frame_data%>%filter(Nonbond<0)
    df_frame_data<-df_frame_data%>%mutate(Elec_interactions="")
    df_frame_data<-df_frame_data%>%mutate(VdW_interactions="")
    df_frame_data$Elec_interactions[df_frame_data$Elec<quantile(df_frame_data$Elec,probs = 0.5)]<-"Elec"
    df_frame_data$VdW_interactions[df_frame_data$VdW<quantile(df_frame_data$VdW,probs = 0.5)]<-"VdW"
    df_frame_data<-df_frame_data%>%filter(Elec<0)
    df_frame_data<-df_frame_data%>%filter(VdW<0)
    df_frame_data<-df_frame_data%>%filter(Nonbond<0)
    df_frame_data<-df_frame_data%>%mutate(interaction=paste0(Elec_interactions,"_",VdW_interactions))
    df_frame_data$interaction[df_frame_data$interaction=="_"]<-"NO"
    df_frame_data$interaction[df_frame_data$Nonbond>0]<-"NO"
    df_frame_data<-df_frame_data%>%group_by(resno,segid)%>%mutate(max=n())
    df_frame_data<-df_frame_data%>%group_by(resno,segid,interaction)%>%mutate(quantity=n())
    df_frame_data<-df_frame_data%>%group_by(resno,segid)%>%mutate(persent=quantity/max)
    df_frame_data<-df_frame_data%>%filter(persent>0.5)
#    ggplot(data=df_frame_data)+
#        geom_point(aes(x=Elec,y=VdW,colour=interaction))+
#        theme_bw()
    df_frame_data<-df_frame_data%>%select(resno, segid, Elec_interactions, 
                                          VdW_interactions,interaction,quantity,max)
    df_frame_data<-unique(df_frame_data)
    df_frame_data<-df_frame_data%>%group_by(resno,segid)%>%mutate(max_quantity=max(quantity))
    df_frame_data<-ungroup(df_frame_data)
    df_frame_data<-df_frame_data%>%filter(quantity==max_quantity)
    
    ggplot(data=df_frame_data)+
        geom_point(aes(x=resno,y=quantity,colour=interaction))+
        theme_bw()
    df_frame_data<-df_frame_data%>%filter(interaction!="NO")
    df_frame_data<-df_frame_data%>%group_by(resno,interaction)%>%mutate(chain_resno=n())
    df_frame_data<-ungroup(df_frame_data)
    df_frame_data<-df_frame_data%>%filter(chain_resno==2)

    ggplot(data=df_frame_data)+
        geom_point(aes(x=resno,y=1,colour=interaction))+
        facet_grid(segid~.)+
        theme_bw()
    df_frame_data<-df_frame_data%>%group_by(resno)%>%mutate(chain_resno=n())
    df_frame_data<-ungroup(df_frame_data)
    df_frame_data<-df_frame_data%>%filter(chain_resno==2)
    df_frame_data<-df_frame_data%>%select(resno,Elec_interactions,VdW_interactions, interaction)
    df_frame_data<-unique(df_frame_data)
    df_frame_data<-df_frame_data%>%group_by(interaction)%>%mutate(amino=paste0(resno,collapse = " "))
    print(df_frame_data$amino)
    df_frame_data<-df_frame_data%>%select(interaction, amino)
    df_frame_data<-unique(df_frame_data)
    
    ggplot(data=df_frame_data)+
        geom_point(aes(x=resno,y=interaction,colour=interaction))+
        #facet_grid(segid~.)+
        theme_bw()
    ggplot(data=df_frame_data)+
        geom_point(aes(x=Elec_interactions,y=VdW_interactions,colour=chain_resno))+
        #facet_grid(segid~.)+
        theme_bw()
    df_frame_data<-df_frame_data%>%group_by(resno,segid)%>%mutate(max=n())
#    df_frame_data<-df_frame_data%>%group_by(resno,segid)%>%mutate(max_interactions=)
    df_frame_data<-ungroup(df_frame_data)

    df_frame_data<-df_frame_data%>%group_by(resno,segid,interaction)%>%mutate(quantity=n())
    df_frame_data<-df_frame_data%>%group_by(resno,segid,interaction)%>%mutate(persent=quantity/max)
    ggplot(data=df_frame_data)+
        geom_point(aes(x=resno,y=persent,colour=interaction))+
        theme_bw()
#    df_frame_data<-df_frame_data%>%filter(Elec_interactions=="Elec")
    df_frame_data<-ungroup(df_frame_data)
    df_frame_data<-df_frame_data%>%select(resno,Elec_interactions,VdW_interactions,interaction)
    df_frame_data<-unique(df_frame_data)
    ggplot(data=df_frame_data)+
        geom_point(aes(x=resno,y=interaction,colour=interaction))+
        theme_bw()
}

