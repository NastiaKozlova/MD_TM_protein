#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
main_part<-c(8)
setwd(part_start)
library(bio3d)
v_parta<-list.files('MD')
v_part<-paste0(part_start,'MD/',v_parta)
v_main<-c('8')
if (!dir.exists(paste0(part_start,'MD_analysis/tcl/'))){dir.create(paste0(part_start,'MD_analysis/tcl/'))}
if (!dir.exists(paste0(part_start,'MD_analysis/hbonds/'))){dir.create(paste0(part_start,'MD_analysis/hbonds/'))}
i<-1
q<-1
j<-1
u<-1
for (j in 1:length(v_parta)) {
    for (q in 1:length(main_part)) {
        part<-paste0(part_start,'MD/',v_parta[j])
        parta<-paste0(v_parta[j])    
        if (!dir.exists(paste0(part,'/din/hbonds_log'))){dir.create(paste0(part,'/din/hbonds_log'))}
        if (!dir.exists(paste0(part,'/din/hbonds_log/',main_part[q],'/'))){dir.create(paste0(part,'/din/hbonds_log/',main_part[q],'/'))}
        if (!dir.exists(paste0(part,'/din/hbonds/'))){dir.create(paste0(part,'/din/hbonds/'))}
        if (!dir.exists(paste0(part,'/din/hbonds/',main_part[q],'/'))){dir.create(paste0(part,'/din/hbonds/',main_part[q],'/'))}
        pdb<-read.pdb(paste0(part,"/din/pdb_second/8/frame_0.pdb"))
        v_segid<-unique(pdb$atom$segid)
        number_frame<-length(list.files(paste0(part,'/din/pdb_second/hbond_',main_part[q],'/')))
        for(u in 1:length(v_segid)){
            df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
            df_tcl[1,1]<-paste('cd', part,"\npackage require hbonds")
            if (number_frame>0){
                for (i in 0:(number_frame-1)) {
                    df_tcl[i+2,1]<-paste('mol new {namd/step5_input.psf} type {psf}')
                    df_tcl[i+2,2]<-paste0('mol addfile {din/pdb_second/hbond_',main_part[q],'/frame_',i,'.pdb} type {pdb}')
                    df_tcl[i+2,3]<-paste0('set protein [atomselect top "protein and segname ',v_segid[u],'" ]')
                    df_tcl[i+2,4]<-paste0('set water [atomselect top "water" ]')
                    df_tcl[i+2,5]<-paste0('hbonds -sel1 $protein -sel2 $water -writefile yes -upsel yes -frames all -dist 3.0 -ang 20 -plot no -outdir din -log hbonds_log/',main_part[q],'/seqid_',v_segid[u],'_frame_',i,'.txt -writefile yes -outfile outfile -polar no -DA both -type all -detailout hbonds/',main_part[q],'/seqid_',v_segid[u],'_frame_',i,'.txt')
                    df_tcl[i+2,6]<-'mol delete all'
                }
                a<-nrow(df_tcl)
                df_tcl[a+1,1]<-'\n\n\nexit now'
                write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta,'_hbond_',main_part[q],"seqid_",v_segid[u],'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
                
                system(command = paste0("vmd -dispdev text -e ",part_start,'MD_analysis/tcl/',parta,'_hbond_',main_part[q],"seqid_",v_segid[u],'.tcl'),ignore.stdout=T,wait = T)
            }
        }
    }
}
#df_fin_conf<-data.frame(matrix(ncol = length(main_part),nrow=length(v_parta)))
#i<-1
#j<-1
#for (i in 1:length(v_parta)) {
#    for (j in 1:length(main_part)) {
#        if (file.exists(paste0(part_start,'MD/',v_parta[i],"/namd/step",v_main[j],".dcd"))){
#            temp_script<-paste0("source ", part_start, "MD_analysis/tcl/",v_parta[i], "_hbond_",      main_part[j],".tcl\n")
#            df_fin_conf[i,j]<-paste0(temp_script)
#        }
#    }
#}
#df_fin_add<-data.frame(matrix(ncol = length(main_part),nrow=1))
#df_fin_add[1,1]<-"\n\nexit now"
#colnames(df_fin_add)<-colnames(df_fin_conf)
#df_fin_conf<-rbind(df_fin_conf,df_fin_add)
#write.table(df_fin_conf,file =paste0(part_start,'MD_analysis/vmd_hbonds_script.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
