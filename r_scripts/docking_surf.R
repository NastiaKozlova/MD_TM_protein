part <- commandArgs(trailingOnly=TRUE)
part_TEMP<-strsplit(part,split = ",",fixed = T)[[1]]
part_start<-part_TEMP[1]
v_search<-part_TEMP[2]
library(ggplot2)
library(bio3d)
library(dplyr)
#setwd(part_start)
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/",v_search,"/")

#system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_group_structure.R ",part_name,",1"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/RMSD_merge_docking_parts.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/calibration_docking_merge_structure.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/merge_docking_parts.R ",part_name,",1"),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/merge_interactions.R ",part_name),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/complex_structure_surf.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/atom_interactions.R ",part_name),ignore.stdout=T,wait = T)
