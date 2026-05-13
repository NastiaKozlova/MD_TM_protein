part_start <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(ggplot2)
library(dplyr)

v_main<-c(8)
main_part<-c(8)
setwd(part_start)

test_10<-seq(from=0,to=1000,by=10)

make_df_ramachadran<-function(part_protein){
    beg<-part_protein$atom$resno[1]
    fin<-part_protein$atom$resno[length(part_protein$atom$resno)]
    tor<-torsion.pdb(part_protein)
    
    df_tor<-data.frame(matrix(nrow = length(tor$phi),ncol=2))
    colnames(df_tor)<-c("phi","psi")
    df_tor[1:length(tor$phi),1]<-tor$phi
    df_tor[1:length(tor$psi),2]<-tor$psi
    df_tor$seq<-NA
    df_tor$seq[1:length(part_protein$atom$resid[atom.select(part_protein,"calpha")$atom])]<-part_protein$atom$resid[atom.select(part_protein,"calpha")$atom]
    df_tor$number<-NA
    df_tor$number<-c(beg:fin)
    df_tor<- df_tor%>% filter(seq!="GLY")
    df_tor<-df_tor%>%mutate(amino=paste(number,seq))
    return(df_tor)
}

filter_df_ramachadran<-function(df_rama,df_rama_start){
    df_rama<-df_rama%>%filter(!is.na(phi))
    df_rama<-df_rama%>%filter(!is.na(psi))
    
    df_rama<-df_rama%>%mutate(favorite=0)
    o<-0
    for (o in 0:(ncol(df_rama_start)/2-1)) {
        
        df_rama_test<-df_rama_start[,(2*o+1):(2*o+2)]
        colnames(df_rama_test)<-c("phi","psi")
        df_rama_test<-df_rama_test%>%filter(!is.na(phi))
        df_rama_test<-df_rama_test%>%filter(!is.na(psi))
        
        for (a in 1:nrow(df_rama)) {
            x0<-df_rama$phi[a]
            y0<-df_rama$psi[a]
            p<-0
            
            for (j in 2:nrow(df_rama_test)) {
                q<- ((df_rama_test$phi[j-1]-x0)*(df_rama_test$phi[j]-x0)+(df_rama_test$psi[j-1]-y0)*(df_rama_test$psi[j]-y0))
                w<- (sqrt((df_rama_test$phi[j-1]-x0)^2+(df_rama_test$psi[j-1]-y0)^2)*sqrt((df_rama_test$phi[j]-x0)^2+(df_rama_test$psi[j]-y0)^2))
                g<-q/w
                if (abs(g)>1) {g<-(1) }
                e<- round(acos(g),digits = 3)
                r<-(-(df_rama_test$phi[j-1]-x0)*(df_rama_test$psi[j]-y0)+(df_rama_test$psi[j-1]-y0)*(df_rama_test$phi[j]-x0))
                
                if (r<0) {p<-p-e  } else{ p<-p+e}
                
                if (abs(p)>df_rama$favorite[a]) {
                    df_rama$favorite[a]<-abs(p)
                }
            }
        } 
    }
    df_rama_filter<-df_rama%>%filter(favorite<6)
    return(df_rama_filter)
}

v_part<-list.files("MD")
if (!dir.exists(paste0(part_start,"MD_analysis/din/"))){dir.create(paste0(part_start,"MD_analysis/din/"))}


for (p in 1:length(v_part)) {
    part<-paste0(part_start,"MD/",v_part[p],"/")
    parta<-paste0(part_start,"MD_analysis/din/",v_part[p],"/")
    if (!dir.exists(paste0(part_start,"MD_analysis/din/",v_part[p]))){dir.create(paste0(part_start,"MD_analysis/din/",v_part[p]))}
    setwd(paste0(part,"din/pdb_second"))
    
    for (main in main_part) {
        if (file.exists(paste0(main,"/frame_",0,".pdb"))){ 
            frame_number<-length(list.files(path = paste0(main)))
            df_topology<-data.frame(matrix(nrow = frame_number, ncol=3))
            colnames(df_topology)<-c("number","frame_number","ramachadran")
            df_topology$number<-0:(frame_number-1)
            df_topology<-df_topology%>%mutate(frame_number=paste0("frame_",number))
            if(!dir.exists(paste0(main,"_rama/"))){(dir.create(paste0(main,"_rama/")))}
            if(!dir.exists(paste0(main,"_prerama/"))){(dir.create(paste0(main,"_prerama/")))}
            for (i in 1:nrow(df_topology)) {
                protein<-read.pdb(paste0(main,"/",df_topology$frame_number[i],".pdb"))
                v_segid<-unique(protein$atom$segid)
                for (q in 1:length(v_segid)) {
                    protein_atom<-atom.select(protein,"protein",segid=v_segid[q])
                    part_protein<-trim.pdb(protein,protein_atom)
                    df_rama<-make_df_ramachadran(part_protein)
                    #df_rama<-df_rama%>%mutate(segid=v_segid[q])
                    df_rama<-df_rama%>%mutate(segid=v_segid[q])
                    #if(!dir.exists(paste0(main,"_rama/")))    { dir.create(paste0(main,"_rama/"))}
                    write.csv(df_rama,file = paste0(main,"_prerama/",df_topology$frame_number[i],"_",v_segid[q],".csv"),row.names = F)
                }
                df_rama<-read.csv(paste0(main,"_prerama/",df_topology$frame_number[i],"_",v_segid[1],".csv"),stringsAsFactors =  F)
                if(length(v_segid)>1){
                    df_rama_add<-read.csv(paste0(main,"_prerama/",df_topology$frame_number[i],"_",v_segid[q],".csv"),stringsAsFactors =  F)
                    df_rama<-rbind(df_rama,df_rama_add)
                }
                write.csv(df_rama,file = paste0(main,"_rama/",df_topology$frame_number[i],".csv"),row.names = F)
            }
            df_rama_start<-read.csv(paste0(part_start,"r_scripts/rama4.csv"))
            i<-1
            q<-1
            for (i in 1:nrow(df_topology)) {
                protein<-read.pdb(paste0(main,"/",df_topology$frame_number[i],".pdb"))
                v_segid<-unique(protein$atom$segid)
                for (q in 1:length(v_segid)) {
                    df_rama<-read.csv(file = paste0(main,"_prerama/",df_topology$frame_number[i],"_",v_segid[q],".csv"),stringsAsFactors =  F)
                    if(!dir.exists(paste0(main,"_ramachad/")))    { dir.create(paste0(main,"_ramachad/"))}
                    if(!dir.exists(paste0(main,"_ramachad_chain/")))    { dir.create(paste0(main,"_ramachad_chain/"))}
                    df_rama_filter<-filter_df_ramachadran(df_rama = df_rama,df_rama_start = df_rama_start)
                    write.csv(df_rama_filter,paste0(main,"_ramachad_chain/",df_topology$frame_number[i],"_",v_segid[q],".csv"),row.names = F)
                }
                df_rama_filter<-read.csv(paste0(main,"_ramachad_chain/",df_topology$frame_number[i],"_",v_segid[1],".csv"),
                                         stringsAsFactors = F)
                df_rama_filter<-df_rama_filter%>%mutate(segid=v_segid[1])
                if(length(v_segid)>1){
                    for (q in 2:length(v_segid)) {
                        df_rama_filter_add<-read.csv(paste0(main,"_ramachad_chain/",df_topology$frame_number[i],"_",v_segid[q],".csv"),
                                                     stringsAsFactors = F) 
                        df_rama_filter_add<-df_rama_filter_add%>%mutate(segid=v_segid[q])
                        df_rama_filter<-rbind(df_rama_filter,df_rama_filter_add)
                    }
                    
                }
                write.csv(df_rama_filter,paste0(main,"_ramachad/",df_topology$frame_number[i],".csv"),row.names = F)
            }
            i<-1
            
            df_rama_filter<-read.csv(paste0(main,"_ramachad/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
            df_rama_filter<-df_rama_filter%>%group_by(segid)%>%summarise(ramachadran=n())
            df_rama_filter<-df_rama_filter%>%mutate(number=df_topology$number[i])%>%
                mutate(frame_number=df_topology$frame_number[i])
            
            for (i in 2:nrow(df_topology)) {
                df_rama_filter_add<-read.csv(paste0(main,"_ramachad/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
                df_rama_filter_add<-df_rama_filter_add%>%group_by(segid)%>%summarise(ramachadran=n())
                df_rama_filter_add<-df_rama_filter_add%>%mutate(number=df_topology$number[i])%>%
                    mutate(frame_number=df_topology$frame_number[i])
                df_rama_filter<-rbind(df_rama_filter,df_rama_filter_add)
            }
            df_rama_filter<-df_rama_filter%>%mutate(time=number/10)
            
            p1<-ggplot(data = df_rama_filter)+
                labs(y="Nimber of aminoasid", x="Time, ns")+
                facet_grid(segid~.)+
                geom_line(aes(x=time,y=ramachadran))+theme_bw()
            ggsave(plot = p1, filename = paste0(parta,main,"_time_Ramachadran.png"),  width = 20, height = 20, units = c("cm"), dpi = 300 )
            write.csv(df_rama_filter,paste0(parta,main,"_time_Ramachadran.csv"),row.names = F)
            df_topology<-df_rama_filter%>%group_by(frame_number)%>%mutate(sum_ramachadran=sum(ramachadran))
            df_topology<-ungroup(df_topology)
            df_topology<-df_topology%>%mutate(ramachadran=sum_ramachadran)
            df_topology<-df_topology%>%group_by(ramachadran)%>%mutate(x_num=n())
            df_topology<-ungroup(df_topology)
            df_topology<-df_topology%>%select(ramachadran, x_num)
            df_topology<-unique(df_topology)
            df_topology<-df_topology%>%mutate(x_num=x_num/frame_number*100)
            write.csv(df_topology,paste0(parta,main,"_Ramachadran.csv"),row.names = F)
        }
        if (file.exists(paste0(parta,main,"_Ramachadran.csv"))){
            df_topology<-read.csv(paste0(parta,main,"_Ramachadran.csv"),stringsAsFactors =  F)
            df_topology<-df_topology%>%arrange(ramachadran)
            df_topology<-df_topology%>%mutate(persent=x_num)
            for (i in 2:nrow(df_topology)) {
                df_topology$persent[i]<-df_topology$persent[i-1]+df_topology$x_num[i]
            }
            df_topology<-df_topology%>%mutate(persent=round(persent,digits = 2))
            df_topology$persent[df_topology$persent>1]<-round(df_topology$persent[df_topology$persent>1],digits = 1)
            p1<-ggplot(data = df_topology)+
                labs(x="Quantity of unfavorite aminoasids",y="persent of time, %")+
                geom_point(aes(x=ramachadran,y=x_num))+
                geom_line(aes(x=ramachadran,y=x_num))+theme_bw()+
                geom_text(aes(x=ramachadran,y=x_num+0.25,label=persent))+
                facet_grid(segid~.)+
                scale_y_continuous(breaks = c(min(round(df_topology$x_num,digits = 0)):max(round(df_topology$x_num,digits = 0))),
                                   labels = c(min(round(df_topology$x_num,digits = 0)):max(round(df_topology$x_num,digits = 0))))+
                scale_x_continuous(breaks = c(min(df_topology$ramachadran):max(df_topology$ramachadran)),
                                   labels = c(min(df_topology$ramachadran):max(df_topology$ramachadran)))
            ggsave(plot = p1, filename = paste0(parta,main,"_Ramachadran.png"),  width = 20, height = 15, units = c("cm"), dpi = 300 )
            df_topology<-read.csv(paste0(parta,main,"_time_Ramachadran.csv"),stringsAsFactors = F)
            df_rama<-read.csv(paste0(main,"_ramachad/",df_topology$frame_number[1],".csv"),stringsAsFactors =  F)
            for (i in 2:nrow(df_topology)) {
                df_rama_filter<-read.csv(paste0(main,"_ramachad/",df_topology$frame_number[i],".csv"),stringsAsFactors = F)
                df_rama<-rbind(df_rama,df_rama_filter)
            }
            df_rama<-df_rama%>%group_by(amino)%>%mutate(x_num=n())
            df_rama<-ungroup(df_rama)
            df_rama<-df_rama%>%select(seq, number ,amino, x_num)
            df_rama<-unique(df_rama) 
            df_rama<-df_rama%>%mutate(x_num=x_num/frame_number*100)
            df_rama_filter<-df_rama
            p1<-ggplot(data = df_rama)+
                labs(x="Nimber of aminoasid", y="persent of time, %")+
                geom_point(aes(x=number,y=x_num))+theme_bw()+
                #geom_line(aes(x=number,y=x_num))+
                facet_grid(segid~.)+
                scale_x_continuous(breaks = test_10,labels = test_10)+
                geom_text(aes(x=number,y=x_num,label=amino),data = df_rama_filter)
            write.csv(df_rama,paste0(parta,main,"_ramachad.csv"),row.names = F)
            ggsave(plot = p1, filename = paste0(parta,main,"_Ramachadran_histogram.png"),  width = 40, height = 20, units = c("cm"), dpi = 300 )
        }
    }
}
