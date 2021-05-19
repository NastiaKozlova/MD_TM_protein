part_start<-'/home/nastia/projects/MD_TM_protein/'
#v_main<-c(6.6,7)
#test disulfid bonds
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepate_test_tcl.R ",part_start),ignore.stdout=T,wait = T) 
#system(command = paste0("vmd -dispdev text -e ",part_start,"vmd_test_script.tcl "),ignore.stdout=T,wait = T) 
#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/disul_bond_leng.R ",part_start),ignore.stdout=T,wait = T)
#errors<-read.csv(paste0(part_start,"Errors_in_disulfid_bonds.csv"),stringsAsFactors = F)
#print(paste0(nrow(errors)," errors in systems\n ",
#             "check disulfid placement bonds in file Errors_in_disulfid_bonds.csv if errors more then 0 and remove systems with wrong disulfid bonds from MD directory"))
#check disulfid placement bonds in file Errors_in_disulfid_bonds.csv and remove systems with wrong disulfid bonds
#MD simulation analysis
if (!dir.exists(paste0(part_start,'MD_analysis/'))){dir.create(paste0(part_start,'MD_analysis/'))}
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/combine_dcd.R ",part_start),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_tcl_din.R ",part_start),ignore.stdout=T,wait = T) 
#vmd
system(command = paste0("vmd -dispdev text -e ",part_start,"MD_analysis/vmd_script.tcl "),ignore.stdout=T,wait = T) 
#check installation of dssp
#sudo apt-get install dssp may help
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/compare_second_str.R ",part_start),ignore.stdout=T,wait = T) 
#rama_4.csv were produse form
#https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/nmr_manual/parameters/manopt_01.html
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/Ramachadran.R ",part_start),ignore.stdout=T,wait = T)

#hbonds
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_XYZ_CA.R ",part_start),ignore.stdout=T,wait = T,show.output.on.console = F)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/prepare_hbonds.R ",part_start),ignore.stdout=T,wait = T)
#vmd
system(command = paste0("vmd -dispdev text -e ",part_start,"vmd_hbonds_script.tcl"),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/test_hbonds.R ",part_start),ignore.stdout=T,wait = T)

#docking
#docking python
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_script.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/convert_pdb_to_pdbqt.R ",part_start),ignore.stdout=T,wait = T)
#for ligands frorm ligand start to ligang
#for receptors from receptor start to receptor
#docking  python
system(command = paste0("chmod +x ",part_start,"docking/script_fin.txt ",part_start),ignore.stdout=T,wait = T)
system(command = paste0(part_start,"docking/script_fin.txt ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/convert_pdbqt_to_pdb.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_pre_analysis.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_group_structure.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking_interactions.R ",part_start),ignore.stdout=T,wait = T)
#ring2 
#check comand ring2 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring2_prepare.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring2_convert.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/ring2_groups.R ",part_start),ignore.stdout=T,wait = T)
#make fin plots
#correct alignemt file for another protein
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/find_conservative_aminoacids.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/make_plots_RMSD_RMSF.R ",part_start),ignore.stdout=T,wait = T)
