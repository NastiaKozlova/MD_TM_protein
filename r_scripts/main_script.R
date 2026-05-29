#readr
part_start<-'path to MD_TM_protein/'
part_start<-paste0(getwd(),"/")
setwd(part_start)
#if you want don't count cout interactions of protein with protein surface surphase_conut<-F
surphase_conut<-T

v_search<-c("center","surf")
#if(!surphase_conut){v_search<-v_search[1]}
v_MD<-list.files(paste0("MD"))
if(!dir.exists("MD_count")){dir.create("MD_count")}
if(!dir.exists("plot")){dir.create("plot")}
i<-1
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_namd_scripts.R ",part_start),ignore.stdout=T,wait = T)

for (i in 1:length(v_MD)) {
  system(command = paste0("cp ",part_start,"MD/",v_MD[i],"/run_namd.txt ",part_start,"MD_count/",v_MD[i],"_run_namd.txt "))
}
#MD simulation analysis
if (!dir.exists(paste0(part_start,'MD_analysis/'))){dir.create(paste0(part_start,'MD_analysis/'))}
#combine all dcd files of Productive MD runs 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/combine_dcd.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_din.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_din_segname.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/Energy_aminoacids_interactions_tcl.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/protein_lipid_interactions.R ",part_start),ignore.stdout=T,wait = T) 

if(file.exists("start/domains_of_interest.csv")){
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_domain_of_interest.R ",part_start),ignore.stdout=T,wait = T) 
}
#check installation of dssp
#sudo apt-get install dssp may help
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/compare_second_str.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/analysis_second_str.R ",part_start),ignore.stdout=T,wait = T)

#rama_4.csv were produse form
#https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/nmr_manual/parameters/manopt_01.html
#install open babel sudo apt-get install openbabel
#download and unzip autodock vina in programs/ http://vina.scripps.edu/download.html
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/Ramachadran.R ",part_start),ignore.stdout=T,wait = T)

#hbonds
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_XYZ_CA.R ",part_start),ignore.stdout=T,wait = T,show.output.on.console = F)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_hbonds.R ",part_start),ignore.stdout=T,wait = T)
#vmd
#system(command = paste0("vmd -dispdev text -e ",part_start,"MD_analysis/vmd_hbonds_script.tcl"),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/test_hbonds.R ",part_start),ignore.stdout=T,wait = T)
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/test_water.R ",part_start),ignore.stdout=T,wait = T)



#docking
#Download programm https://ccsb.scripps.edu/mgltools/downloads/
#docking python
#/home/nastia/projects/MD_TM_protein/r_scripts/docking/docking_prepare_receptor_pdb.R

#make fin plots
#correct alignemt file for another protein
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_conservative_aminoacids.R ",part_start),ignore.stdout=T,wait = T)
#collect MD simulation data in minimal amount of dataframes 

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plots_RMSD_RMSF_frame.R ",part_start),ignore.stdout=T,wait = T)
#rewrite
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plots_RMSD_RMSF.R ",part_start),ignore.stdout=T,wait = T)

#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/claster_analysis_frame_data.R ",part_start),ignore.stdout=T,wait = T)


system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_prepare_receptor_pdb.R ",part_start),ignore.stdout=T,wait = T)
v_search<-c("center","surf")
i<-1
#part_name<-paste0(part_start,",",v_search[i])
for (i in 1:length(v_search)) {
    part_name<-paste0(part_start,",",v_search[i])
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/copy_files_for_docking.R ",part_name),ignore.stdout=T,wait = T)
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/r_scripts/docking_add_serf_active_center.R ",part_start),ignore.stdout=T,wait = T)
}
for (i in 1:length(v_search)) {   
    part_name<-paste0(part_start,",",v_search[i])
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_premain.R ",part_name),ignore.stdout=T,wait = T)
}
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_add_serf_active_center.R ",part_start),ignore.stdout=T,wait = T)
#for (i in 1:length(v_list_proteins)) {
    for (j in 1:length(v_search)) {
        part<-paste0(part_start,",",v_search[j])
        system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_script.R ",part),ignore.stdout=T,wait = T)
        
        print("Run")
        print(paste0(part,"script/script_readme.txt"))
    }
#}
j<-2
for (j in 1:length(v_search)) {
    part<-paste0(part_start,",",v_search[j])
    #part<-paste0(part_start,"MD_analysis/docking/docking_first/",v_search[j],"/")
    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_main.R ",part),ignore.stdout=T,wait = T)
} 
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/",v_search,"/")
for (j in 1:length(v_search)) {

    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_group_structure.R ",part_name[j],",",1),ignore.stdout=T,wait = T)
}
#if you want don't count cout interactions of protein with protein serfuce v_surphase_conut<-F  
part_name<-paste0(part_start,",",v_search)
for (j in 1:length(v_search)) {

    system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_surf.R ",part_name[j]),ignore.stdout=T,wait = T)
}
#if(surphase_conut){
#  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_surf.R ",part_start),ignore.stdout=T,wait = T)
for (j in 1:length(v_search)) {
    part_name<-paste0(part_start,"MD_analysis/docking/docking_first/",v_search,"/")
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_convert_active_center.R ",part_start),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_active_center.R ",part_start),ignore.stdout=T,wait = T)
  system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/cascade_docking.R ",part_start),ignore.stdout=T,wait = T)
}

#ring2 
#check comand ring2 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring/ring2_prepare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring/ring2_convert.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring/ring2_groups.R ",part_start),ignore.stdout=T,wait = T)
#ligands_plasement
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/ligands_plasement.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/ligand_placement_analysis.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_ligand_path.R ",part_start),ignore.stdout=T,wait = T)

#ligand_placement_analysis
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_plot.R ",part_start),ignore.stdout=T,wait = T)
arrange_ligand_positions
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/ligand_placement.R ",part_start),ignore.stdout=T,wait = T)

# make random plot
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/count comparition.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/domain_interactions.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/claster_analysis_frame_data.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/claster_analysis_frame_data_RMSD.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/claster_analysis_frame_data_all.R ",part_start),ignore.stdout=T,wait = T)

#test_clusterisation
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/test_clusterisation.R ",part_start),ignore.stdout=T,wait = T)


# system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/test_clusterisation.R ",part_start),ignore.stdout=T,wait = T)

#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/frame_statistical_analysis_two_systems.R ",part_start),ignore.stdout=T,wait = T)
