part_start = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(bio3d)
library(dplyr)
setwd(part_start)
v_parta<-list.files("ligand_start")


for (j in 1:length(v_parta)) {
  b<-strsplit(v_parta[j],split = ".",fixed = T)[[1]]
  pdb_name<-b
  system(command = paste0("obabel ",part_start,"ligand_start/",pdb_name, ".pdb -O ",part_start,"ligand/",pdb_name, ".pdbqt"),ignore.stdout=T,wait = T)
}

