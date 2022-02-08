part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(dplyr)
library(bio3d)
library(GGally)
library(network)
test_10<-seq(from=0,to=1000,by=10)
df_topology<-read.csv(paste0(part_start,"start/df_topology.csv"),stringsAsFactors =  F)
part<-paste0(part_start,"MD_analysis/")
setwd(part)
parta<-paste0(part_start,"MD/")
df_all_systems<-read.csv(paste0(part_start,"start/all_systems.csv"),stringsAsFactors = F)

df_all_systems$system_name<-as.character(df_all_systems$system_name)
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
main_part<-c(8)
i<-nrow(df_all_systems)
main<-main_part[2]
p<-1

q<-2

if(!dir.exists(paste0(part,"fin_data/TMD_interactions/"))){dir.create(paste0(part,"fin_data/TMD_interactions/"))}
if(!dir.exists(paste0(part,"fin_data/domain_importance/"))){dir.create(paste0(part,"fin_data/domain_importance/"))}
#if(!dir.exists(paste0(part,"fin_data/aminoacids_interactions/"))){dir.create(paste0(part,"fin_data/aminoacids_interactions/"))}
#if(!dir.exists(paste0(part,"fin_data/aminoacids_importance/"))){dir.create(paste0(part,"fin_data/aminoacids_importance/"))}
if(!dir.exists(paste0(part,"fin_plots/aminoacids_interactions/"))){dir.create(paste0(part,"fin_plots/aminoacids_interactions/"))}
if(!dir.exists(paste0(part,"fin_plots/aminoacids_interactions_sort/"))){dir.create(paste0(part,"fin_plots/aminoacids_interactions_sort/"))}
if(!dir.exists(paste0(part,"fin_plots/TMD_interactions/"))){dir.create(paste0(part,"fin_plots/TMD_interactions/"))}


mm.col <- c("ECD"= "#ff69b4",
            "TMD"= "#4783B7",
            "CD" = "#F48037", 
            "loop"= "#80B682")
q<-1
for (q in 1:nrow(df_all_systems)) {
  if(file.exists(paste0("fin_data/str_data/",df_all_systems$fin_name[q],".csv"))){
    if(file.exists(paste0("din/",df_all_systems$fin_name[q],"/8_ring_interaction.csv"))){
      print(df_all_systems$fin_name[q])
      df_seq<-read.csv(paste0("fin_data/str_data/",df_all_systems$fin_name[q],".csv"),stringsAsFactors = F)
      df_ring<-read.csv(paste0("fin_data/aminoacids_interactions/",df_all_systems$fin_name[q],".csv"),stringsAsFactors = F)
      df_ring<-df_ring%>%filter(persent_intractions>0)
      df_ring<-left_join(df_ring,df_seq,by=c("number.x"="resno"))
      df_ring<-left_join(df_ring,df_seq,by=c("number.y"="resno"))

      df_ring_sorted<-df_ring%>%filter(persent_intractions>quantile(df_ring$persent_intractions,probs = 0.75))
      df_ring_sorted<-df_ring_sorted%>%filter(abs(number.y-number.x)>5)
      df_ring<-df_ring%>%filter(abs(number.y-number.x)>5)
      df_interactions_TMD<-df_ring_sorted%>%filter(number.x<number.y)
      df_interactions_TMD<-df_interactions_TMD%>%select( type.x,topology.x, type.y,topology.y)
      df_interactions_TMD<-df_interactions_TMD%>%mutate(TMD_bond=paste(type.x,type.y))%>%group_by(TMD_bond)%>%mutate(number_TMD_intractions=n())
      df_interactions_TMD<-ungroup(df_interactions_TMD)
      df_interactions_TMD<-df_interactions_TMD%>%filter(type.x!=type.y)
      df_interactions_TMD<-unique(df_interactions_TMD)
#      write.csv(df_interactions_TMD,paste0("fin_data/TMD_interactions/",df_all_systems$fin_name[q],".csv"),
#                row.names = F)
      df_interactions_TMD<-df_interactions_TMD%>%select(type.x,type.y,number_TMD_intractions)
      df_TMD_importance<-df_ring_sorted%>%select(type.x,topology.x)
      df_TMD_importance<-df_TMD_importance%>%group_by(type.x)%>%mutate(number_TMD_intractions=n())
#      write.csv(df_TMD_importance,paste0("fin_data/domain_importance/",df_all_systems$fin_name[q],".csv"),
#                row.names = F)

      if(nrow(df_interactions_TMD)>1){
        df_TMD_importance<-unique(df_TMD_importance)
        rownames(df_TMD_importance) <- df_TMD_importance$type.x
        vlist<-list(edges =df_interactions_TMD,vertices =df_TMD_importance)
        
       # rownames(vlist$vertices) <- vlist$vertices$type.x
        mm.net <- network(vlist$edges[, 1:3], directed = T)
        mm.net %v% "topology.x" <- as.character(
          vlist$vertices[ network.vertex.names(mm.net), "topology.x"]
        )
        # type color palette
        # create plot for ggnet2
        set.seed(10052006)
        p<-ggnet2(mm.net,#color = mm.col[ mm.net %v% "topology.x" ],
                  label = T,
                  edge.size = "number_TMD_intractions",
                  size = 20, mode = "kamadakawai", label.size = 5)

        ggsave(p,filename = paste0("fin_plots/TMD_interactions/",df_all_systems$fin_name[q],"_",df_all_systems$Membrane[q],"_",df_all_systems$Structure[q],".png"), width = 30, height = 30, units = c("cm"), dpi = 200 ) 
      }
      if(nrow(df_ring)>1){

        df_ring<-df_ring%>%filter(number.x<number.y)
        rownames(df_seq) <- df_seq$resno
        vlist<-list(edges =df_ring,vertices =df_seq)
        
#        rownames(vlist$vertices) <- vlist$vertices$resno
        mm.net <- network(vlist$edges[, 1:3], directed = FALSE)
        mm.net %v% "topology" <- as.character(
          vlist$vertices[ network.vertex.names(mm.net), "topology"]
        )

        p<-ggnet2(mm.net, color = mm.col[ mm.net %v% "topology" ],
                  label = T,
                  edge.size = "persent_intractions",
                  size = 5, vjust = -0.6, mode = "fruchtermanreingold", label.size = 5)
        ggsave(p,filename = paste0("fin_plots/aminoacids_interactions/",df_all_systems$fin_name[q],"_",df_all_systems$Membrane[q],"_",df_all_systems$Structure[q],".png"), width = 30, height = 30, units = c("cm"), dpi = 200 ) 
      }
      if(nrow(df_ring_sorted)>1){

        df_ring_sorted<-df_ring_sorted%>%filter(number.x<number.y)
        rownames(df_ring_sorted) <- df_ring_sorted$resno
        vlist<-list(edges =df_ring_sorted,vertices =df_seq)
        
#        rownames(vlist$vertices) <- vlist$vertices$resno
        mm.net <- network(vlist$edges[, 1:3], directed = FALSE)
        mm.net %v% "topology" <- as.character(
          vlist$vertices[ network.vertex.names(mm.net), "topology"]
        )
        p<-ggnet2(mm.net, color = mm.col[ mm.net %v% "topology" ],
                  label = T,
                  edge.size = "persent_intractions",
                  node.size = "importance_ring", vjust = -0.6, mode = "fruchtermanreingold", label.size = 5)
        ggsave(p,filename = paste0("fin_plots/aminoacids_interactions_sort/",df_all_systems$fin_name[q],"_",df_all_systems$Membrane[q],"_",df_all_systems$Structure[q],".png"), width = 30, height = 30, units = c("cm"), dpi = 200 ) 
      }
    }
  }
}  
