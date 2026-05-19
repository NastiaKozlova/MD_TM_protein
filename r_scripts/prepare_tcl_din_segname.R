#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
#install.packages("stringr")
library(stringr)
library(bio3d)
library(dplyr)
setwd(part_start)
v_parta<-list.files('MD')
v_part<-paste0(part_start,'MD/',v_parta)
namd_exe<-paste0(part_start,"programs/NAMD_2.14_Linux-x86_64-multicore/namd2 ")

if (!dir.exists(paste0(part_start,'MD_analysis/tcl/'))){dir.create(paste0(part_start,'MD_analysis/tcl/'))}
if (!dir.exists(paste0(part_start,'MD_analysis/hbonds/'))){dir.create(paste0(part_start,'MD_analysis/hbonds/'))}
i<-1
q<-1
j<-1
p<-1
for (j in 1:length(v_parta)) {
    part<-paste0(part_start,'MD/',v_parta[j],"/")
    if(file.exists(paste0(part,"namd/step",8,".dcd"))){
        parta<-paste0(v_parta)
        
        if (!dir.exists(paste0(part,'/din'))){dir.create(paste0(part,'/din'))}
        if (!dir.exists(paste0(part,'/din/pdb_second'))){dir.create(paste0(part,'/din/pdb_second'))}
        if (!dir.exists(paste0(part,'/din/pdb_second/',8))){dir.create(paste0(part,'/din/pdb_second/',8))}
        if (!dir.exists(paste0(part,'/din/pdb_second/hbond_',8))){dir.create(paste0(part,'/din/pdb_second/hbond_',8))}
        
        if (!dir.exists(paste0(part,'/din/RMSD'))){dir.create(paste0(part,'/din/RMSD'))}
        if (!dir.exists(paste0(part,'/din/RMSF'))){dir.create(paste0(part,'/din/RMSF'))}
        if (!dir.exists(paste0(part,'/din/SASA'))){dir.create(paste0(part,'/din/SASA'))}
        if (!dir.exists(paste0(part,'/din/Energy'))){dir.create(paste0(part,'/din/Energy'))}
        if (!dir.exists(paste0(part,'/din/Energy_protein_lipid'))){dir.create(paste0(part,'/din/Energy_protein_lipid'))}
        if (!dir.exists(paste0(part,'/din/hbonds'))){dir.create(paste0(part,'/din/hbonds'))}
        
        df_conf<-read.table(file = paste0(part,"/namd/step7_production.inp"),sep="?")
        v_paraneters<-c()
        for (i in 1:nrow(df_conf)) {
            if(length(df_conf[i,1])>0&(grepl( df_conf[i,1], pattern = "parameters", fixed = TRUE))){
                v_paraneters<-c(v_paraneters,df_conf[i,1])
            }
        }
        a<-  str_locate(v_paraneters, "toppar/")
        i<-1
        for (i in 1:length(v_paraneters)) {
            seq<-strsplit(v_paraneters[i],split = "",fixed = T)[[1]]
            v_paraneters[i]<-paste0(seq[a[i,1]:length(seq)],collapse = "")
        }
        v_paraneters<-paste0("-par namd/",v_paraneters)
        v_paraneters<-paste(v_paraneters,collapse = " ")
        
        pdb<-read.pdb(paste0(part,"din/pdb_second/8/frame_0.pdb"))
        v_segid<-pdb$atom$segid
        for (p in 1:length(v_segid)){
            df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
            df_tcl[1,1]<-paste('cd', part,"\npackage require namdenergy")
            df_tcl[2,1]<-paste0('mol new {namd/step5_input.psf} type {psf}')
            df_tcl[3,1]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
            df_tcl[4,1]<-paste0('set sel0 [atomselect top "water"]')
            df_tcl[5,1]<-paste0('set sel1 [atomselect top "lipid or resname POPG"]')
            df_tcl[6,1]<-paste0('set sel2 [atomselect top "protein and segname ',v_segid[p],'"]')
            
            df_tcl[7,1]<-paste0('namdenergy -sel $sel2  $sel0 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/protein_water_energy_segname_',v_segid[p],'.txt -switch 10 -exe ',namd_exe,' ',v_paraneters)
            df_tcl[8,1]<-paste0('namdenergy -sel $sel2  $sel1 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/protein_lipid_energy_segname_',v_segid[p],'.txt -switch 10 -exe ',namd_exe,' ',v_paraneters)
            df_tcl[9,1]<-paste0('namdenergy -sel $sel2        -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/protein__segname_',v_segid[p],'.txt -switch 10 -exe ',namd_exe,' ',v_paraneters)
            df_tcl[10,1]<-'mol delete all\n\n\n exit now'
            write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_Energy_',v_segid[p],'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
            
            system(command = paste0("vmd -dispdev text -e ",part_start,'MD_analysis/tcl/',parta[j],'_Energy_',v_segid[p],'.tcl'),ignore.stdout=T,wait = T)
        }
        
    }
}