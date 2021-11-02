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
df_all_systems<-df_all_systems%>%mutate(Progress=="notdone")
df_all_systems<-df_all_systems%>%mutate(fin_name=paste0("charmm-gui-",system_name))
v_part<-list.files(paste0(part,"din"))
main_part<-c(8)
i<-nrow(df_all_systems)
main<-main_part[2]
p<-1

q<-2
#v_topology<-c()
#for (i in 1:nrow(df_topology)) {
#  a<-c(df_topology$seq_beg[i]:df_topology$seq_end[i]) 
#  v_topology<-c(v_topology,a)
#}
#v_resno<-df_pdb$resno
#v_resno<-v_resno[!v_resno%in%v_topology]
if(!dir.exists(paste0(part,"fin_data/TMD_interactions/"))){dir.create(paste0(part,"fin_data/TMD_interactions/"))}
if(!dir.exists(paste0(part,"fin_data/domain_importance/"))){dir.create(paste0(part,"fin_data/domain_importance/"))}
if(!dir.exists(paste0(part,"fin_data/aminoacids_interactions/"))){dir.create(paste0(part,"fin_data/aminoacids_interactions/"))}
if(!dir.exists(paste0(part,"fin_data/aminoacids_importance/"))){dir.create(paste0(part,"fin_data/aminoacids_importance/"))}
if(!dir.exists(paste0(part,"fin_plots/aminoacids_interactions/"))){dir.create(paste0(part,"fin_plots/aminoacids_interactions/"))}
if(!dir.exists(paste0(part,"fin_plots/TMD_interactions/"))){dir.create(paste0(part,"fin_plots/TMD_interactions/"))}


mm.col <- c("ECD"= "#ff69b4",
            "TMD"= "#4783B7",
            "CD" = "#F48037", 
            "loop"= "#80B682")

for (q in 1:nrow(df_all_systems)) {
  if(file.exists(paste0("docking/receptor_start/",df_all_systems$fin_name[q],".pdb"))){
    if(file.exists(paste0("din/",df_all_systems$fin_name[q],"/8_ring_interaction.csv"))){
      df_ring<-read.csv(paste0("din/",df_all_systems$fin_name[q],"/8_ring_interaction.csv"),stringsAsFactors = F)
      df_ramachadran<-read.csv(paste0("din/",df_all_systems$fin_name[q],"/8_time_Ramachadran.csv"),stringsAsFactors = F)
      pdb<-read.pdb(paste0("docking/receptor_start/",df_all_systems$fin_name[q],".pdb"))
      df_pdb<-pdb$atom
      df_pdb<-df_pdb%>%filter(elety=="CA")
      df_pdb<-df_pdb%>%mutate(topology=NA)
      for (i in 1:nrow(df_topology)) {
        df_pdb$type[(df_pdb$resno>=df_topology$seq_beg[i])&(df_pdb$resno<=df_topology$seq_end[i])]<-df_topology$type[i]
        df_pdb$topology[(df_pdb$resno>=df_topology$seq_beg[i])&(df_pdb$resno<=df_topology$seq_end[i])]<-df_topology$topology[i]
      }
      df_pdb<-df_pdb%>%filter(type!="ATOM")
      df_ring<-df_ring%>%filter((number.y-number.x)>5)
      df_ring<-df_ring%>%filter(all>0)
      df_ring<-df_ring%>%filter(all>nrow(df_ramachadran)*0.8)
      df_ring<-df_ring[df_ring$number.x%in%df_pdb$resno,]
      df_ring<-df_ring[df_ring$number.y%in%df_pdb$resno,]
      df_ring1<-df_ring%>%select(number.x,number.y)
      df_ring2<-df_ring%>%select(number.y,number.x)
      colnames(df_ring2)<-colnames(df_ring1)
      df_ring<-rbind(df_ring1,df_ring2)
      df_interactions<-df_ring1
      df_interactions<-unique(df_interactions)
      df_interactions<-left_join(df_interactions,df_pdb,by=c("number.x"="resno"))
      df_interactions<-left_join(df_interactions,df_pdb,by=c("number.y"="resno"))
      df_pdb_amino<-df_pdb%>%select(resno, type,topology)
      df_pdb<-df_pdb%>%select( type,topology)
      df_pdb<-unique(df_pdb)
      df_interactions<-df_interactions%>%select(type.x, type.y)

      df_interactions<-df_interactions%>%filter(type.x!=type.y)
      df_interaction<-df_interactions%>%mutate(bond=paste0(type.x,type.y))
      df_interactions<-df_interaction%>%group_by(bond)%>%mutate(number=n())
      df_interactions<-ungroup(df_interactions)
      df_interactions<-df_interactions%>%select(type.x,type.y,number)
      df_interactions<-unique(df_interactions)
      
      df_interaction1<-df_interaction
      df_interaction2<-df_interaction%>%select(type.y,type.x,bond)
      colnames(df_interaction2)<-colnames(df_interaction1)
      df_interaction<-rbind(df_interaction1,df_interaction2)
      df_interaction<-df_interaction%>%group_by(type.x)%>%mutate(number=n())
      df_interaction<-df_interaction%>%select(type.x,number)
      df_interaction<-unique(df_interaction)
      df_interaction<-ungroup(df_interaction)
      colnames(df_interaction)<-c("type", "number")
      df_pdb<-left_join(df_pdb,df_interaction,by="type")
      
      write.csv(df_interactions,paste0("fin_data/TMD_interactions/",df_all_systems$fin_name[q],".csv"),
                row.names = F)
      write.csv(df_pdb,paste0("fin_data/domain_importance/",df_all_systems$fin_name[q],".csv"),
                row.names = F)
      df_pdb<-df_pdb%>%mutate(number=number/2)
      df_pdb<-unique(df_pdb)

      df_interaction_amino<-df_ring%>%group_by(number.x)%>%mutate(importance=n())
      df_interaction_amino<-ungroup(df_interaction_amino)
      df_pdb_amino<-left_join(df_pdb_amino,df_interaction_amino,by=c("resno"="number.x"))
      #df_pdb_amino<-df_pdb_amino%>%filter(!is.na(importance))
      df_pdb_amino$number.y<-NULL
      df_pdb_amino<-unique(df_pdb_amino)
      
      df_interactions_bonds<-df_ring%>%select("number.x","number.y")
      df_interactions_bonds<-df_interactions_bonds%>%filter(number.x<number.y)
      

      write.csv(df_interactions_bonds,paste0("fin_data/aminoacids_interactions/",df_all_systems$fin_name[q],".csv"),
                row.names = F)
      df_pdb_amino$importance[is.na(df_pdb_amino$importance)]<-0
      df_pdb_amino<-df_pdb_amino%>%arrange(desc(importance))
      write.csv(df_pdb_amino,paste0("fin_data/aminoacids_importance/",df_all_systems$fin_name[q],".csv"),
                row.names = F)
      if(nrow(df_interactions)>1){
        vlist<-list(edges =df_interactions,vertices =df_pdb)
        
        rownames(vlist$vertices) <- vlist$vertices$type
        mm.net <- network(vlist$edges[, 1:3], directed = FALSE)
        mm.net %v% "topology" <- as.character(
          vlist$vertices[ network.vertex.names(mm.net), "topology"]
        )
        # type color palette

        # create plot for ggnet2
        set.seed(10052006)
        p<-ggnet2(mm.net, color = mm.col[ mm.net %v% "topology" ],
                  label = T,
                  edge.size = "number",
                  size = 20, mode = "kamadakawai", label.size = 5)
        ggsave(p,filename = paste0("fin_plots/TMD_interactions/",df_all_systems$fin_name[q],"_",df_all_systems$disulfid_bonds[q],"_",df_all_systems$Glyco[q],".png"), width = 30, height = 30, units = c("cm"), dpi = 200 ) 

        vlist<-list(edges =df_interactions_bonds,vertices =df_pdb_amino)
        
        rownames(vlist$vertices) <- vlist$vertices$resno
        mm.net <- network(vlist$edges[, 1:2], directed = FALSE)
        mm.net %v% "topology" <- as.character(
          vlist$vertices[ network.vertex.names(mm.net), "topology"]
        )
        
        p<-ggnet2(mm.net, color = mm.col[ mm.net %v% "topology" ],
               label = T,
               #label.color = mm.col[ mm.net %v% "topology" ],
               size = 5, vjust = -0.6, mode = "fruchtermanreingold", label.size = 5)
        ggsave(p,filename = paste0("fin_plots/aminoacids_interactions/",df_all_systems$fin_name[q],"_",df_all_systems$disulfid_bonds[q],"_",df_all_systems$Glyco[q],".png"), width = 30, height = 30, units = c("cm"), dpi = 200 ) 
        
      }
    }
  }
}  
