#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
#install.packages("stringr")
library(stringr)
setwd(part_start)
df_domains_of_interest<-read.csv("start/domains_of_interest.csv",stringsAsFactors = F)
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
    if (!dir.exists(paste0(part,'/din/hbonds'))){dir.create(paste0(part,'/din/hbonds'))}
    for (p in 1:nrow(df_domains_of_interest)) {
      
      
      
      df_conf<-read.table(file = paste0(part,"/namd/step7_production.inp"),sep="?")
      v_paraneters<-c()
      for (i in 1:nrow(df_conf)) {
        if(length(df_conf[i,1])>0&(grepl( df_conf[i,1], pattern = "parameters", fixed = TRUE))){
          v_paraneters<-c(v_paraneters,df_conf[i,1])
        }
      }
      #str_view(v_paraneters, "toppar/")
      a<-  str_locate(v_paraneters, "toppar/")
      i<-1
      for (i in 1:length(v_paraneters)) {
        seq<-strsplit(v_paraneters[i],split = "",fixed = T)[[1]]
        v_paraneters[i]<-paste0(seq[a[i,1]:length(seq)],collapse = "")
      }
      v_paraneters<-paste0("-par namd/",v_paraneters)
      v_paraneters<-paste(v_paraneters,collapse = " ")
      
      
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
      df_tcl[1,1]<-paste('cd', part,"\npackage require namdenergy")
      df_tcl[2,1]<-paste0('mol new {namd/step5_input.psf} type {psf}')
      df_tcl[3,1]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
      df_tcl[4,1]<-paste0('set sel0 [atomselect top "water"]')
      df_tcl[5,1]<-paste0('set sel1 [atomselect top "lipid or resname POPG or resname POPS"]')
      df_tcl[6,1]<-paste0('set sel2 [atomselect top "protein and ',
                          'resid>',(df_domains_of_interest$start[p]-1),
                          ' and resid<',(df_domains_of_interest$finish[p]+1),'"]')
      
      df_tcl[7,1]<-paste0('namdenergy -sel $sel2  $sel0 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/',df_domains_of_interest$domain_name[p],'_water_energy_',8,'.txt -switch 10 -exe ',namd_exe,' ',v_paraneters)
      df_tcl[8,1]<-paste0('namdenergy -sel $sel2  $sel1 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/',df_domains_of_interest$domain_name[p],'_lipid_energy_',8,'.txt -switch 10 -exe ',namd_exe,' ',v_paraneters)
      df_tcl[9,1]<-paste0('namdenergy -sel $sel2        -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/',df_domains_of_interest$domain_name[p],'_',8,             '.txt -switch 10 -exe ',namd_exe,' ',v_paraneters)
      df_tcl[10,1]<-'mol delete all\n\n\n exit now'
      write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],"_",df_domains_of_interest$domain_name[p],'_Energy_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
      system(command = paste0("vmd -dispdev text -e ",part_start,'MD_analysis/tcl/',parta[j],"_",df_domains_of_interest$domain_name[p],'_Energy_',8,'.tcl'),ignore.stdout=T,wait = T)
      
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 13))
      df_tcl[1,1]<-paste('cd', part,'\nmol new {namd/step5_input.psf} type {psf}')
      df_tcl[1,2]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
      df_tcl[1,3]<-paste0('set protein_0 [atomselect top "protein and name CA  and ',
                          'resid>',(df_domains_of_interest$start[p]-1),
                          ' and resid<',(df_domains_of_interest$finish[p]+1),'" frame 0]')
      df_tcl[1,4]<-paste0('set n [molinfo top get numframes]')
      df_tcl[1,5]<-paste0('set output [open din/RMSD/',8,"_",df_domains_of_interest$domain_name[p],'.txt w] ')
      df_tcl[1,6]<-paste0('for {set i 0} {$i < $n} {incr i} {')
      df_tcl[1,7]<-paste0('set protein_i [atomselect top "protein and name CA and ',
                          'resid>',(df_domains_of_interest$start[p]-1),
                          ' and resid<',(df_domains_of_interest$finish[p]+1),'" frame $i]')
      df_tcl[1,8]<-paste0('set rmsd [measure rmsd $protein_0 $protein_i]')
      df_tcl[1,9]<-paste0('puts $output "$i $rmsd"')
      df_tcl[1,10]<-paste0('}')
      df_tcl[1,11]<-paste0('puts "output file: $n din/RMSD/',8,"_",df_domains_of_interest$domain_name[p],'.txt"')
      df_tcl[1,12]<-paste0('close $output')
      df_tcl[1,13]<-paste0('mol delete all\n\n\n exit now')
      write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],"_",df_domains_of_interest$domain_name[p],'_RMSD_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
      system(command = paste0("vmd -dispdev text -e ",part_start,'MD_analysis/tcl/',parta[j],"_",df_domains_of_interest$domain_name[p],'_RMSD_',8,'.tcl'),ignore.stdout=T,wait = T)
      
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
      df_tcl[1,1]<-paste('cd', part)
      df_tcl[2,1]<-paste0('mol new {namd/step5_input.psf} type {psf}')
      df_tcl[3,1]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
      df_tcl[4,1]<-paste0('set protein [atomselect top "protein and ',
                          'resid>',(df_domains_of_interest$start[p]-1),
                          ' and resid<',(df_domains_of_interest$finish[p]+1),'"]')
      df_tcl[5,1]<-paste0('set n [molinfo top get numframes]')
      df_tcl[6,1]<-paste0('set output [open din/SASA/',8,"_",df_domains_of_interest$domain_name[p],'.txt w] ')
      df_tcl[7,1]<-paste0('for {set i 0} {$i < $n} {incr i} {')
      df_tcl[8,1]<-paste0('molinfo top set frame $i')
      df_tcl[9,1]<-paste0('set sasa_protein [measure sasa 1.4 $protein]')
      df_tcl[10,1]<-paste0('puts $output " $i $sasa_protein "')
      df_tcl[11,1]<-paste0('}')
      df_tcl[12,1]<-paste0('puts "output file: $n din/SASA/',8,"_",df_domains_of_interest$domain_name[p],'.txt"')
      df_tcl[13,1]<-paste0('close $output')
      df_tcl[14,1]<-paste0('mol delete all\n\n\n exit now')
      write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],"_",df_domains_of_interest$domain_name[p],'_SASA_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
      system(command = paste0("vmd -dispdev text -e ",part_start,'MD_analysis/tcl/',parta[j],"_",df_domains_of_interest$domain_name[p],'_SASA_',8,'.tcl'),ignore.stdout=T,wait = T)
    }
  }
}
