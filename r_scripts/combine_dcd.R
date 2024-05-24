#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
v_parta<-list.files('MD')
v_part<-paste0(part_start,'MD/',v_parta)

if (!dir.exists(paste0(part_start,'MD_analysis/tcl/'))){dir.create(paste0(part_start,'MD_analysis/tcl/'))}
if (!dir.exists(paste0(part_start,'MD_analysis/hbonds/'))){dir.create(paste0(part_start,'MD_analysis/hbonds/'))}
i<-1
q<-1
j<-1
num_model<-1000
for (j in 1:length(v_parta)) {
  
  part<-paste0(part_start,'MD/',v_parta[j])
  parta<-paste0(v_parta[j])
  
  if (!dir.exists(paste0(part,'/din'))){dir.create(paste0(part,'/din'))}
  if (!dir.exists(paste0(part,'/din/pdb_second'))){dir.create(paste0(part,'/din/pdb_second'))}
  
  if (!dir.exists(paste0(part,'/din/RMSD'))){dir.create(paste0(part,'/din/RMSD'))}
  if (!dir.exists(paste0(part,'/din/RMSF'))){dir.create(paste0(part,'/din/RMSF'))}
  if (!dir.exists(paste0(part,'/din/SASA'))){dir.create(paste0(part,'/din/SASA'))}
  if (!dir.exists(paste0(part,'/din/Energy'))){dir.create(paste0(part,'/din/Energy'))}
  if (!dir.exists(paste0(part,'/din/hbonds'))){dir.create(paste0(part,'/din/hbonds'))}
  
  #combine and align MD simulation dcd files
  df_tcl<-data.frame(matrix(nrow = 2,ncol = 2))
  df_tcl[1,1]<-paste0('cd ', part,'/namd/\n\npackage require animate\n')
  df_tcl[1,2]<-paste0('mol new {step5_input.psf} type {psf}')
  for (p in 1:num_model) {
    if(file.exists(paste0(part,'/namd/step7.',p,'_production.dcd'))){
      df_tcl[p+1,1]<-paste0('mol addfile {',part,'/namd/step7.',p,'_production.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
    }
  }
  df_tcl[p+2,1]<-paste0('set ref [atomselect top "protein and name CA" frame 0]\n',
                        'set sel [atomselect top "protein and name CA"]\n',
                        'set all [atomselect top all]\n',
                        'set n [molinfo top get numframes]\n',
                        'set fin [expr $n-1]\n',
                        'for { set i 1 } { $i < $n } { incr i } {\n',
                        '  $sel frame $i\n',
                        '  $all frame $i\n',
                        '  $all move [measure fit $sel $ref]\n',
                        '}\n')
  df_tcl[p+3,1]<-paste0('animate write dcd step8.dcd waitfor all')
  df_tcl[p+4,1]<-paste0('mol delete all\n\nexit now')
  
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta,'_combine.tcl'),sep = ' ',na = '' ,row.names = F,col.names = F,quote = F)
  print( paste0('vmd -dispdev text -e ',part_start,'MD_analysis/tcl/',parta,'_combine.tcl'))
  if(file.exists(paste0(part,'/namd/step7.1_production.dcd'))){
    system(command = paste0('vmd -dispdev text -e ',part_start,'MD_analysis/tcl/',parta,'_combine.tcl'),ignore.stdout=T,wait = T) 
  }
  #combine MD simulation dcd files
  df_tcl<-data.frame(matrix(nrow = 2,ncol = 2))
  df_tcl[1,1]<-paste0('cd ', part,'/namd/\n\npackage require animate\n')
  df_tcl[1,2]<-paste0('mol new {step5_input.psf} type {psf}')
  for (p in 1:num_model) {
    if(file.exists(paste0(part,'/namd/step7.',p,'_production.dcd'))){
      df_tcl[p+1,1]<-paste0('mol addfile {',part,'/namd/step7.',p,'_production.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
    }
  }
  df_tcl[p+2,1]<-paste0('animate write dcd step8_nonalign.dcd waitfor all')
  df_tcl[p+3,1]<-paste0('mol delete all\n\nexit now')
  
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta,'_combine_non_align.tcl'),sep = ' ',na = '' ,row.names = F,col.names = F,quote = F)
  print( paste0('vmd -dispdev text -e ',part_start,'MD_analysis/tcl/',parta,'_combine_non_align.tcl'))
  if(file.exists(paste0(part,'/namd/step7.1_production.dcd'))){
    system(command = paste0('vmd -dispdev text -e ',part_start,'MD_analysis/tcl/',parta,'_combine_non_align.tcl'),ignore.stdout=T,wait = T) 
  }
  df_tcl<-data.frame(matrix(nrow = 2,ncol = 2))
  df_tcl[1,1]<-paste0('cd ', part,'/namd/\n\npackage require animate\n')
  df_tcl[1,2]<-paste0('mol new {step5_input.psf} type {psf}\n')
  df_tcl[1,3]<-paste0('mol addfile {',part,'/namd/step8.dcd} type {dcd} first 0 last -1 step 1 waitfor all\n')
  df_tcl[1,4]<-paste0('set s1 [atomselect top "protein"]\n\n')
  df_tcl[1,5]<-paste0('animate write dcd {',part,'/namd/step8_protein.dcd} waitfor all sel $s1\n\n')  
  df_tcl[1,6]<-paste0('[atomselect top "protein" frame 0] writepdb ',part,'/namd/PCA_pdb.pdb\n\n')
  df_tcl[1,7]<-paste0('mol delete all\n\nexit now')
  
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta,'_protein_align.tcl'),sep = ' ',na = '' ,row.names = F,col.names = F,quote = F)
  print( paste0('vmd -dispdev text -e ',part_start,'MD_analysis/tcl/',parta,'_protein_align.tcl'))
  if(file.exists(paste0(part,'/namd/step7.1_production.dcd'))){
    system(command = paste0('vmd -dispdev text -e ',part_start,'MD_analysis/tcl/',parta,'_protein_align.tcl'),ignore.stdout=T,wait = T) 
  }
}