#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
v_parta<-list.files('MD')
v_part<-paste0(part_start,'MD/',v_parta)
namd_exe<-paste0(part_start,"programs/NAMD_2.14_Linux-x86_64-multicore/namd2 ")

if (!dir.exists(paste0(part_start,'MD_analysis/tcl/'))){dir.create(paste0(part_start,'MD_analysis/tcl/'))}
if (!dir.exists(paste0(part_start,'MD_analysis/hbonds/'))){dir.create(paste0(part_start,'MD_analysis/hbonds/'))}
i<-1
q<-1
for (j in 1:length(v_parta)) {
  part<-paste0(part_start,'MD/',v_parta[j])
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
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste('cd', part,"\npackage require namdenergy")
  df_tcl[2,1]<-paste0('mol new {namd/step5_input.psf} type {psf}')
  df_tcl[3,1]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[4,1]<-paste0('set sel0 [atomselect top "water"]')
  df_tcl[5,1]<-paste0('set sel1 [atomselect top "lipid or resname POPG"]')
  df_tcl[6,1]<-paste0('set sel2 [atomselect top "protein"]')
  df_tcl[7,1]<-paste0('namdenergy -sel $sel2  $sel0 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/protein_water_energy_',8,'.txt -switch 10 -exe ',namd_exe,' -par namd/toppar/par_all36_cgenff.prm -par namd/toppar/toppar_all36_na_rna_modified.str -par namd/toppar/toppar_all36_label_fluorophore.str -par namd/toppar/par_all36m_prot.prm -par namd/toppar/par_all36_carb.prm -par namd/toppar/toppar_all36_carb_imlab.str -par namd/toppar/toppar_all36_synthetic_polymer.str -par namd/toppar/toppar_all36_label_spin.str -par namd/toppar/par_interface.prm -par namd/toppar/toppar_all36_prot_c36m_d_aminoacids.str -par namd/toppar/par_all36_na.prm -par namd/toppar/toppar_all36_prot_modify_res.str -par namd/toppar/toppar_all36_lipid_cholesterol.str -par namd/toppar/toppar_all36_lipid_sphingo.str -par namd/toppar/par_all36_lipid.prm -par namd/toppar/toppar_all36_carb_glycolipid.str -par namd/toppar/toppar_all36_lipid_lps.str -par namd/toppar/toppar_all36_prot_retinol.str -par namd/toppar/toppar_all36_carb_glycopeptide.str -par namd/toppar/toppar_all36_prot_na_combined.str -par namd/toppar/toppar_all36_prot_heme.str -par namd/toppar/toppar_all36_na_nad_ppi.str -par namd/toppar/toppar_all36_lipid_ether.str -par namd/toppar/toppar_water_ions.str -par namd/toppar/toppar_all36_prot_arg0.str -par namd/toppar/toppar_ions_won.str -par namd/toppar/toppar_all36_lipid_bacterial.str -par namd/toppar/toppar_all36_nano_lig.str -par namd/toppar/toppar_all36_lipid_detergent.str -par namd/toppar/toppar_all36_prot_fluoro_alkanes.str -par namd/toppar/toppar_all36_lipid_prot.str -par namd/toppar/toppar_all36_lipid_miscellaneous.str -par namd/toppar/toppar_all36_lipid_hmmm.str -par namd/toppar/toppar_all36_lipid_yeast.str -par namd/toppar/toppar_dum_noble_gases.str')
  df_tcl[8,1]<-paste0('namdenergy -sel $sel2  $sel1 -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/protein_lipid_energy_',8,'.txt -switch 10 -exe ',namd_exe,' -par namd/toppar/par_all36_cgenff.prm -par namd/toppar/toppar_all36_na_rna_modified.str -par namd/toppar/toppar_all36_label_fluorophore.str -par namd/toppar/par_all36m_prot.prm -par namd/toppar/par_all36_carb.prm -par namd/toppar/toppar_all36_carb_imlab.str -par namd/toppar/toppar_all36_synthetic_polymer.str -par namd/toppar/toppar_all36_label_spin.str -par namd/toppar/par_interface.prm -par namd/toppar/toppar_all36_prot_c36m_d_aminoacids.str -par namd/toppar/par_all36_na.prm -par namd/toppar/toppar_all36_prot_modify_res.str -par namd/toppar/toppar_all36_lipid_cholesterol.str -par namd/toppar/toppar_all36_lipid_sphingo.str -par namd/toppar/par_all36_lipid.prm -par namd/toppar/toppar_all36_carb_glycolipid.str -par namd/toppar/toppar_all36_lipid_lps.str -par namd/toppar/toppar_all36_prot_retinol.str -par namd/toppar/toppar_all36_carb_glycopeptide.str -par namd/toppar/toppar_all36_prot_na_combined.str -par namd/toppar/toppar_all36_prot_heme.str -par namd/toppar/toppar_all36_na_nad_ppi.str -par namd/toppar/toppar_all36_lipid_ether.str -par namd/toppar/toppar_water_ions.str -par namd/toppar/toppar_all36_prot_arg0.str -par namd/toppar/toppar_ions_won.str -par namd/toppar/toppar_all36_lipid_bacterial.str -par namd/toppar/toppar_all36_nano_lig.str -par namd/toppar/toppar_all36_lipid_detergent.str -par namd/toppar/toppar_all36_prot_fluoro_alkanes.str -par namd/toppar/toppar_all36_lipid_prot.str -par namd/toppar/toppar_all36_lipid_miscellaneous.str -par namd/toppar/toppar_all36_lipid_hmmm.str -par namd/toppar/toppar_all36_lipid_yeast.str -par namd/toppar/toppar_dum_noble_gases.str')
  df_tcl[9,1]<-paste0('namdenergy -sel $sel2        -vdw -elec -nonb -cutoff 12 -skip 0 -ofile din/Energy/protein_',8,             '.txt -switch 10 -exe ',namd_exe,' -par namd/toppar/par_all36_cgenff.prm -par namd/toppar/toppar_all36_na_rna_modified.str -par namd/toppar/toppar_all36_label_fluorophore.str -par namd/toppar/par_all36m_prot.prm -par namd/toppar/par_all36_carb.prm -par namd/toppar/toppar_all36_carb_imlab.str -par namd/toppar/toppar_all36_synthetic_polymer.str -par namd/toppar/toppar_all36_label_spin.str -par namd/toppar/par_interface.prm -par namd/toppar/toppar_all36_prot_c36m_d_aminoacids.str -par namd/toppar/par_all36_na.prm -par namd/toppar/toppar_all36_prot_modify_res.str -par namd/toppar/toppar_all36_lipid_cholesterol.str -par namd/toppar/toppar_all36_lipid_sphingo.str -par namd/toppar/par_all36_lipid.prm -par namd/toppar/toppar_all36_carb_glycolipid.str -par namd/toppar/toppar_all36_lipid_lps.str -par namd/toppar/toppar_all36_prot_retinol.str -par namd/toppar/toppar_all36_carb_glycopeptide.str -par namd/toppar/toppar_all36_prot_na_combined.str -par namd/toppar/toppar_all36_prot_heme.str -par namd/toppar/toppar_all36_na_nad_ppi.str -par namd/toppar/toppar_all36_lipid_ether.str -par namd/toppar/toppar_water_ions.str -par namd/toppar/toppar_all36_prot_arg0.str -par namd/toppar/toppar_ions_won.str -par namd/toppar/toppar_all36_lipid_bacterial.str -par namd/toppar/toppar_all36_nano_lig.str -par namd/toppar/toppar_all36_lipid_detergent.str -par namd/toppar/toppar_all36_prot_fluoro_alkanes.str -par namd/toppar/toppar_all36_lipid_prot.str -par namd/toppar/toppar_all36_lipid_miscellaneous.str -par namd/toppar/toppar_all36_lipid_hmmm.str -par namd/toppar/toppar_all36_lipid_yeast.str -par namd/toppar/toppar_dum_noble_gases.str')
  df_tcl[10,1]<-'mol delete all'
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_Energy_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
  
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 13))
  df_tcl[1,1]<-paste('cd', part,'\nmol new {namd/step5_input.psf} type {psf}')
  df_tcl[1,2]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[1,3]<-paste0('set protein_0 [atomselect top "protein and name CA" frame 0]')
  df_tcl[1,4]<-paste0('set n [molinfo top get numframes]')
  df_tcl[1,5]<-paste0('set output [open din/RMSD/',8,'.txt w] ')
  df_tcl[1,6]<-paste0('for {set i 0} {$i < $n} {incr i} {')
  df_tcl[1,7]<-paste0('set protein_i [atomselect top "protein and name CA" frame $i]')
  df_tcl[1,8]<-paste0('set rmsd [measure rmsd $protein_0 $protein_i]')
  df_tcl[1,9]<-paste0('puts $output "$i $rmsd"')
  df_tcl[1,10]<-paste0('}')
  df_tcl[1,11]<-paste0('puts "output file: $n din/RMSD/',8,'.txt"')
  df_tcl[1,12]<-paste0('close $output')
  df_tcl[1,13]<-paste0('mol delete all')
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_RMSD_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 12))
  df_tcl[1,1]<-paste('cd', part,'\nmol new {namd/step5_input.psf} type {psf}')
  df_tcl[1,2]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[1,3]<-paste0('set protein [atomselect top "protein and name CA"]')
  df_tcl[1,4]<-paste0('set n [molinfo top get numframes]')
  df_tcl[1,5]<-paste0('set output [open din/RMSF/',8,'.txt w] ')
  df_tcl[1,6]<-paste0('set rmsf [measure rmsf $protein]')
  df_tcl[1,7]<-paste0('foreach x $rmsf {')
  df_tcl[1,8]<-paste0('puts $output $x')
  df_tcl[1,9]<-paste0('}')
  df_tcl[1,10]<-paste0('puts "output file: $n din/RMSF/',8,'.txt"')
  df_tcl[1,11]<-paste0('close $output')
  df_tcl[1,12]<-paste0('mol delete all')
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_RMSF_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
  
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste('cd', part)
  df_tcl[2,1]<-paste0('mol new {namd/step5_input.psf} type {psf}')
  df_tcl[3,1]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[4,1]<-paste0('set protein [atomselect top "protein"]')
  df_tcl[5,1]<-paste0('set n [molinfo top get numframes]')
  df_tcl[6,1]<-paste0('set output [open din/SASA/',8,'.txt w] ')
  df_tcl[7,1]<-paste0('for {set i 0} {$i < $n} {incr i} {')
  df_tcl[8,1]<-paste0('molinfo top set frame $i')
  df_tcl[9,1]<-paste0('set sasa_protein [measure sasa 1.4 $protein]')
  df_tcl[10,1]<-paste0('puts $output " $i $sasa_protein "')
  df_tcl[11,1]<-paste0('}')
  df_tcl[12,1]<-paste0('puts "output file: $n din/SASA/',8,'.txt"')
  df_tcl[13,1]<-paste0('close $output')
  df_tcl[14,1]<-paste0('mol delete all')
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_SASA_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
  df_tcl[1,1]<-paste('cd', part,'\nmol new {namd/step5_input.psf} type {psf}')
  df_tcl[1,2]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[1,3]<-paste0('set nf [molinfo top get numframes]')
  df_tcl[1,4]<-paste0('for {set i 0 } {$i < $nf} {incr i} {')
  df_tcl[1,5]<-paste0('[atomselect top "protein" frame $i] writepdb din/pdb_second/',8,'/frame_$i.pdb')
  df_tcl[1,6]<-paste0('}')
  df_tcl[1,7]<-'mol delete all'
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_Second_str_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
  
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
  df_tcl[1,1]<-paste('cd', part,'\nmol new {namd/step5_input.psf} type {psf}')
  df_tcl[1,2]<-paste0('mol addfile {namd/step',8,'.dcd} type {dcd} first 0 last -1 step 1 waitfor all')
  df_tcl[1,3]<-paste0('set nf [molinfo top get numframes]')
  df_tcl[1,4]<-paste0('for {set i 0 } {$i < $nf} {incr i} {')
  df_tcl[1,5]<-paste0('[atomselect top all frame $i] writepdb din/pdb_second/hbond_',8,'/frame_$i.pdb')
  df_tcl[1,6]<-paste0('}')
  df_tcl[1,7]<-'mol delete all'
  write.table(df_tcl,file =paste0(part_start,'MD_analysis/tcl/',parta[j],'_Second_str_hbond_',8,'.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)
}

df_fin_conf<-data.frame(matrix(ncol = length(8),nrow=length(v_parta)))
i<-1
j<-1
for (i in 1:length(v_parta)) {
  for (j in 1:length(8)) {
    if (file.exists(paste0(part_start,'MD/',v_parta[i],"/namd/step",8,".dcd"))){
      temp_script<-paste0("source ", part_start, "MD_analysis/tcl/",  v_parta[i], "_Energy_",     8,".tcl\n",
                          "source ", part_start, "MD_analysis/tcl/",  v_parta[i], "_RMSD_",       8, ".tcl\n",
                          "source ", part_start, "MD_analysis/tcl/",  v_parta[i], "_RMSF_",       8, ".tcl\n",
                          "source ", part_start, "MD_analysis/tcl/",  v_parta[i], "_SASA_",       8, ".tcl\n",
                          "source ", part_start, "MD_analysis/tcl/",  v_parta[i], "_Second_str_",        8, ".tcl\n",
                          "source ", part_start, "MD_analysis/tcl/",   v_parta[i], "_Second_str_hbond_", 8,".tcl\n")
      df_fin_conf[i,j]<-paste0(temp_script)
    }
  }
}
df_fin_add<-data.frame(matrix(ncol = length(8),nrow=1))
df_fin_add[1,1]<-"exit now"
colnames(df_fin_add)<-colnames(df_fin_conf)
df_fin_conf<-rbind(df_fin_conf,df_fin_add)
write.table(df_fin_conf,file =paste0(part_start,'MD_analysis/vmd_script.tcl'),sep = '\n', quote = F,na = '' ,row.names = F,col.names = F)

