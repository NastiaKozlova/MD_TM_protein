part_start <- commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)

part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/docking_group_structure.R ",part_name,",1"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/RMSD_merge_docking_center.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/calibration_merge_structure_center.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/merge_docking_center.R ",part_name,",1"),ignore.stdout=T,wait = T)


system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/merge_interactions_center.R ",part_name,",1"),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/complex_structure_center.R ",part_name),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/docking/atom_interactions_center.R ",part_name),ignore.stdout=T,wait = T)