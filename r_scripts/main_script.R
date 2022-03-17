#readr
part_start<-'path to MD_TM_protein/'
part_start<-paste0(getwd(),"/")
setwd(part_start)

v_MD<-list.files(paste0("MD"))
if(!dir.exists("MD_count")){dir.create("MD_count")}
if(!dir.exists("plot")){dir.create("plot")}
i<-1
for (i in 1:length(v_MD)) {
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_namd_scripts.R ",part_start,"MD/",v_MD[i],"/"),ignore.stdout=T,wait = T)
  system(command = paste0("cp ",part_start,"MD/",v_MD[i],"/run_namd.txt ",part_start,"MD_count/",v_MD[i],"_run_namd.txt "))
}
#MD simulation analysis
if (!dir.exists(paste0(part_start,'MD_analysis/'))){dir.create(paste0(part_start,'MD_analysis/'))}
#combine all dcd files of Productive MD runs 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/combine_dcd.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_din.R ",part_start),ignore.stdout=T,wait = T) 
if(file.exists("start/domains_of_interest.csv")){
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_domain_of_interest.R ",part_start),ignore.stdout=T,wait = T) 
}
#check installation of dssp
#sudo apt-get install dssp may help
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/compare_second_str.R ",part_start),ignore.stdout=T,wait = T) 
#rama_4.csv were produse form
#https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/nmr_manual/parameters/manopt_01.html
#install open babel sudo apt-get install openbabel
#download and unzip autodock vina in programs/ http://vina.scripps.edu/download.html
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/Ramachadran.R ",part_start),ignore.stdout=T,wait = T)

#hbonds
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_XYZ_CA.R ",part_start),ignore.stdout=T,wait = T,show.output.on.console = F)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_hbonds.R ",part_start),ignore.stdout=T,wait = T)
#vmd
system(command = paste0("vmd -dispdev text -e ",part_start,"MD_analysis/vmd_hbonds_script.tcl"),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/test_hbonds.R ",part_start),ignore.stdout=T,wait = T)

#docking
#Download programm https://ccsb.scripps.edu/mgltools/downloads/
#docking python
#/home/nastia/projects/MD_TM_protein/r_scripts/docking/docking_prepare_receptor_pdb.R
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_prepare_receptor_pdb.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_add_serf_active_center.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_main.R ",part_start),ignore.stdout=T,wait = T)

#ring2 
#check comand ring2 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring/ring2_prepare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring/ring2_convert.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring/ring2_groups.R ",part_start),ignore.stdout=T,wait = T)
#make fin plots
#correct alignemt file for another protein
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_conservative_aminoacids.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plots_RMSD_RMSF.R ",part_start),ignore.stdout=T,wait = T)
# make random plot
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/count comparition.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/domain_interactions.R ",part_start),ignore.stdout=T,wait = T)
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/amino_analysis.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/fin_structure_prepare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/atom_interactions.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/claster_analysis_frame_data.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/frame_statistical_analysis_two_systems.R ",part_start),ignore.stdout=T,wait = T)
