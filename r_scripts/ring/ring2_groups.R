part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
library(readr)
test_10<-seq(from=0,to=1000,by=10)

setwd(part_start)
v_part<-list.files("MD")
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/ring2/"))){dir.create(paste0(part_start,"MD_analysis/ring2/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/ring2/script"))){dir.create(paste0(part_start,"MD_analysis/ring2/script"))}
p<-1
i<-1
y<-1
main_part<-c(8)
a<-list.files(paste0(part_start,"start/sequense/"))
seq<-read.fasta(paste0(part_start,"start/sequense/",a))

main<-main_part[2]
p<-30
for (p in 1:length(v_part)) {
  part<-paste0(part_start,"MD/",v_part[p],"/din/pdb_second/")
  parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
  setwd(part)
  for (main in main_part) {
    if(file.exists(paste0(parta,main,"_time_Ramachadran.csv"))){
      df_topology<-read.csv(paste0(parta,main,"_time_Ramachadran.csv"),stringsAsFactors = F)
      df_seq<-t(seq$ali)
      colnames(df_seq)<-"seq"
      df_seq<-as.data.frame(df_seq)
      df_seq<-df_seq%>%mutate(number=1:nrow(df_seq))
      df_seq<-df_seq%>%mutate(C="C")
      df_seq_bond<-full_join(df_seq,df_seq,by="C")
      df_seq_bond<-df_seq_bond%>%filter(number.x<number.y)
      df_seq_bond<-df_seq_bond%>%mutate(bond=paste(number.x,number.y,sep = "-"))
      df_seq_bond<-df_seq_bond%>%mutate(all=0)
      df_seq_bond<-df_seq_bond%>%mutate(HBOND=0)
      df_seq_bond<-df_seq_bond%>%mutate(VDW=0)
      df_seq_bond<-df_seq_bond%>%mutate(PIPISTACK=0)
      df_seq_bond<-df_seq_bond%>%mutate(IONIC=0)
      df_seq_bond<-df_seq_bond%>%mutate(IAC=0)
      df_seq_bond<-df_seq_bond%>%mutate(PICATION=0)
      df_seq_bond$C<-NULL
      
      df_topology<-df_topology%>%select(number,frame_number)
      df_topology<-df_topology%>%mutate(all=NA)
      df_topology<-df_topology%>%mutate(HBOND=NA)
      df_topology<-df_topology%>%mutate(VDW=NA)
      df_topology<-df_topology%>%mutate(PIPISTACK=NA)
      df_topology<-df_topology%>%mutate(IONIC=NA)
      df_topology<-df_topology%>%mutate(IAC=NA)
      df_topology<-df_topology%>%mutate(PICATION=NA)
#      rm(df_seq)
      for (i in 1:nrow(df_topology)) {
        df_ring<-read.csv(paste0(main,"_ring_mod/",df_topology$frame_number[i],".txt"),stringsAsFactors = F)  
        df_ring_HBOND<-df_ring%>%filter(Interaction=="HBOND")
        df_topology$HBOND[i]<-nrow(df_ring_HBOND) 
      
        df_ring_VDW<-df_ring%>%filter(Interaction=="VDW")
        df_topology$VDW[i]<-nrow(df_ring_VDW) 
      
        df_ring_PIPISTACK<-df_ring%>%filter(Interaction=="PIPISTACK")
        df_topology$PIPISTACK[i]<-nrow(df_ring_PIPISTACK) 
      
        df_ring_IONIC<-df_ring%>%filter(Interaction=="IONIC")
        df_topology$IONIC[i]<-nrow(df_ring_IONIC) 
      
        df_ring_IAC<-df_ring%>%filter(Interaction=="IAC")
        df_topology$IAC[i]<-nrow(df_ring_IAC) 
      
        df_ring_PICATION<-df_ring%>%filter(Interaction=="PICATION")
        df_ring<-    df_ring%>%select(NodeId1, NodeId2, bond )
        df_ring<-    unique(df_ring)
        df_topology$all[i]<-nrow(df_ring) 
      
        df_topology$PICATION[i]<-nrow(df_ring_PICATION) 
      
        for (j in 1:nrow(df_ring)) {
          df_seq_bond$all[df_seq_bond$bond==df_ring$bond[j]]<-df_seq_bond$all[df_seq_bond$bond==df_ring$bond[j]]+1
        }
        if (nrow(df_ring_HBOND)>0) {
          for (j in 1:nrow(df_ring_HBOND)) {
            df_seq_bond$HBOND[df_seq_bond$bond==df_ring_HBOND$bond[j]]<-df_seq_bond$HBOND[df_seq_bond$bond==df_ring_HBOND$bond[j]]+1
          }
        }
        if (nrow(df_ring_VDW)>0) {
          for (j in 1:nrow(df_ring_VDW)) {
            df_seq_bond$VDW[df_seq_bond$bond==df_ring_VDW$bond[j]]<-df_seq_bond$VDW[df_seq_bond$bond==df_ring_VDW$bond[j]]+1
          }
        }
        if (nrow(df_ring_PIPISTACK)>0) {
          for (j in 1:nrow(df_ring_PIPISTACK)) {
            df_seq_bond$PIPISTACK[df_seq_bond$bond==df_ring_PIPISTACK$bond[j]]<-df_seq_bond$PIPISTACK[df_seq_bond$bond==df_ring_PIPISTACK$bond[j]]+1
          }
        }
        if (nrow(df_ring_IONIC)>0) {
          for (j in 1:nrow(df_ring_IONIC)) {
            df_seq_bond$IONIC[df_seq_bond$bond==df_ring_IONIC$bond[j]]<-df_seq_bond$IONIC[df_seq_bond$bond==df_ring_IONIC$bond[j]]+1
          }
        }
        if (nrow(df_ring_IAC)>0) {
          for (j in 1:nrow(df_ring_IAC)) {
            df_seq_bond$IAC[df_seq_bond$bond==df_ring_IAC$bond[j]]<-df_seq_bond$IAC[df_seq_bond$bond==df_ring_IAC$bond[j]]+1
          }
        }
        if (nrow(df_ring_PICATION)>0) {
          for (j in 1:nrow(df_ring_PICATION)) {
            df_seq_bond$PICATION[df_seq_bond$bond==df_ring_PICATION$bond[j]]<-df_seq_bond$PICATION[df_seq_bond$bond==df_ring_PICATION$bond[j]]+1
          }
        }
      }
      print(nrow(df_seq_bond))
      write.csv(df_seq_bond,paste0(parta,main,"_ring_interaction.csv"),row.names = F)
    }
  }
}
