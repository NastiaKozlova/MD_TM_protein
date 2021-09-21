# MD_TM_protein
Group of scripts with could help run and analyse long trajectories of MD simulations protein-membrane systems _**only by NAMD**_ prepared by CHARMM GUI

CHARMM GUI output should be placed into _**MD_TM_protein/MD**_ folder
Scripts are counting protein parameters such as SASA, Ramachandran, protein energy and could interact amino acids of protein with water. 
They also can trace the network of interprotein interactions. Scripts also could analyse could ligands interact with ligands. 

This group of scripts supports preparing and analyzing MD simulations of the complex receptor(globular proteins)-ligand.
Scripts are separate to the four parts:

1. Protein structure stabilization by MD simulation using NAMD
2. Receptor-ligand complex constraction

The main script from with program will run is _*r_scripts/master_script.R*_

## Nessesary files

list of membrane structures in pdb into _*MD_TMD_protein/start/all_systems.csv*_
Example of the table

| system\_name | Progress       | Membrane | Structure |
| ------------ | -------------- | -------- | --------- |
| 2322184819   | to do          | POPG     | WT        |
| 2322253581   | to do          | POPE     | WT        |

Sequence of protein into _*MD_TMD_protein/start/sequence/Protein_name.fasta*_
Sequences of proteins to align into _*MD_TMD_protein/start/alignment/Protein_name.fasta*_
Protein topology into _*MD_TMD_protein/start/df_topology.csv*_
List of active center aminoacids _*MD_TMD_protein/start/active_center.csv*_
Example of the table:
| type   | amino | resno |
| ------ | ----- | ----- |
| center | PHE   | 27    |
| center | TRP   | 151   |
Type distinguish different centers of interest
Ligand structures in pdb into _*MD_TMD_protein/start/ligand_start*_




# Protein structure stablilisation using MD simulations (NAMD)

CHARMM GUI output should be placed into _**MD_TM_protein/MD**_ folder
" MD_TM_protein/r_scripts/make_namd_scripts.R" -- prepare all scripts and structures to run MD simulations. 
MD simulations are contained equvlibration according to CHARMM-GUI manual and  1000 ns productive simulation(1000\*1ns).
You can adjast length of simulation by editing _*num_din*_ parameter. _*num_din*_ time in ns

After this step you can find your namd scripts to _*MD_TMD_protein/MD_count/*_
This script has to be altered to accommodate different configurations of computers and run them OUTSIDE of the R environment.


# Prediction of structure receptor-ligand complex this docking AutoDock

You should put all of ligand structures in pdb into _*MD_TMD_protein/start/sequence/lacY.fasta*_
You should put all of ligand structures in pdb into _*MD_TMD_protein/start/active_center.csv*_
You should put all of ligand structures in pdb into _*MD_TMD_protein/start/ligand_start*_
You should put all of ligand structures in pdb into _*MD_TMD_protein/start/alignment/lacY.fasta*_
You should put all of ligand structures in pdb into _*MD_TMD_protein/start/df_topology.csv*_
You should put all of ligand structures in pdb into _*MD_TMD_protein/start/all_systems.csv*_

# Receptor-ligand structure stablilisation using MD simulations (NAMD)

The main script from with program will run is _*r_scripts/master_script.R*_
"r_scripts/prepare_to_stabilisation_MD.R" -- prepare all scripts and structures 
to run MD simulations. MD simulations are contained 0.1ns minimisation, 0.3ns heating,
2 ns eqvilibration and 100 ns productive simulation(100\*1ns). 
You can alter time by changing _*num_din*_

After this step you can find your namd script _*MD\_globular\_protein/r\_scripts/namd\_script.txt*_
This script has to be altered to accommodate different configurations of computers. 


# Dependencise

## Programs

1. NAMD (NAMD could be downliaded from https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=NAMD you should copy namd multicore into _**MD_TM_protein/programm/**_ folder)
2. VMD (VMD could be downliaded https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD)
3. R (sudo apt-get update sudo apt-get install r-base r-base-dev or read https://cran.r-project.org/bin/linux/debian/)
4. dssp Debian and Ubuntu (sudo apt install dssp)
5. ring2 (could be downloaded form http://old.protein.bio.unipd.it/download/)
6. Autodock (Autodock vina could be downliaded _**from http://vina.scripps.edu/download.html **_ you should copy Autodock into _**MD_TM_protein/programm/autodock_vina_1_1_2_linux_x86**_ folder)
7. MGLTools (MGLTools into could be downliaded from https://ccsb.scripps.edu/mgltools/downloads/, and you shold chouse this version of the program _*mgltools_Linux-x86_64_1.5.7_Install (Linux 64 GUI installer 109Mb)*_ or rename folder atherthords to _* MD_TM_protein/programs/MGLTools-1.5.7*_)

## R packages

1. bio3d install.packages("bio3d")
2. cowlpot install.packages("cowlpot")
3. ggplot2 install.packages("ggplot2")
4. dplyr install.packages("dplyr")
5. httr install.packages("httr")

# References 

1. J. Lee, D.S. Patel, J. Ståhle, S-J. Park, N.R. Kern, S. Kim, J. Lee, X. Cheng, M.A. Valvano, O. Holst, Y. Knirel, Y. Qi, S. Jo, J.B. Klauda, G. Widmalm, and W. Im (2019) CHARMM-GUI Membrane Builder for Complex Biological Membrane Simulations with Glycolipids and Lipoglycans. J. Chem. Theory Comput. 15:775-786 https://pubs.acs.org/doi/10.1021/acs.jctc.8b01066
2. NAMD 
3. VMD Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, vol. 14, pp. 33-38. doi:10.1016/0263-7855(96)00018-5 http://linkinghub.elsevier.com/retrieve/pii/0263785596000185
4. R
5. dssp 
 A)A series of PDB related databases for everyday needs. Wouter G Touw, Coos Baakman, Jon Black, Tim AH te Beek, E Krieger, Robbie P Joosten, Gert Vriend. Nucleic Acids Research 2015 January; 43(Database issue): D364-D368. 10.1093/nar/gku1028 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4383885/

2) Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Kabsch W, Sander C, Biopolymers. 1983 22 2577-2637. PMID: 6667333; UI: 84128824. https://pubmed.ncbi.nlm.nih.gov/6667333/ doi: 10.1002/bip.360221211
6. ring2 Piovesan D., Minervini G., Tosatto S.C.E. The RING 2.0 web server for high quality residue interaction networks. (2016) Nucleic Acids Research, 44 (W1), pp. W367-74.
7. Autodock  http://doi.wiley.com/10.1002/jcc.21256
8. MGLTools
9. bio3d
10. cowlpot
11. ggplot2
12. dplyr
13. httr
