part_start <- commandArgs(trailingOnly=TRUE)
library(dplyr)
library(ggplot2)
library(bio3d)
part_name<-list.files(paste0(part_start,"MD_analysis/"))

#v_seq<-read.csv(paste0(part_start,"start/sequence/",seq_name[1]))

part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")
#/media/akozlova/Data2/Nastia_projects/lacY/MD_TMD_protein/MD_analysis/docking/docking_first/
setwd(part_name)
df_merge<-read.csv(file = paste0(part_name,"din/df_structure.csv"),stringsAsFactors = F)
df_convert<-df_merge%>%select( receptor,system_name,Membrane,Structure)
df_convert<-unique(df_convert)

if (!dir.exists(paste0(part_name,"tost/"))) {dir.create(paste0(part_name,"tost/"))}
if (!dir.exists(paste0(part_name,"plot/"))) {dir.create(paste0(part_name,"plot/"))}
if (!dir.exists(paste0(part_name,"all_interactions/"))) {dir.create(paste0(part_name,"all_interactions/"))}
if (!dir.exists(paste0(part_start,"MD_analysis/fin_data/docking_data/"))) {
  dir.create(paste0(part_start,"MD_analysis/fin_data/docking_data/"))}
if (!dir.exists(paste0(part_name,"din/tost/"))) {
  dir.create(paste0(part_name,"din/tost/"))}
i<-1
for (i in 1:nrow(df_merge)) {
  if (!dir.exists(paste0(part_name,"din/tost/",df_merge$receptor[i],"_",df_merge$ligand[i]))) {
    dir.create(paste0(part_name,"din/tost/",df_merge$receptor[i],"_",df_merge$ligand[i]))}
  receptor_name<-paste0(part_name,"receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_name,"din/str_fin/",df_merge$name.x[i])
  pdb_receptor<-read.pdb(receptor_name)
  pdb_ligand<-read.pdb(ligand_name)
  pdb_complex<-cat.pdb(pdb_receptor, pdb_ligand, rechain=TRUE)
  
  #  pdb<-read.pdb(paste0(part_name,"str_fin/",df_merge$name.x[i]))
  write.pdb(pdb_complex,paste0(part_name,"din/tost/",df_merge$receptor[i],"_",df_merge$ligand[i],"/",df_merge$structure_order[i],"_",df_merge$name.x[i]))
}
df_select<-df_merge%>%select(receptor,ligand)
df_select<-unique(df_select)
df_merge<-df_merge%>%mutate(folder_name="surf")
for (i in 1:nrow(df_merge)) {
  if(strsplit(df_merge$name.x[i],split = "_")[[1]][3]=="serf"){
    df_merge$folder_name[i]<-"interaction_serf"
  }else{
    df_merge$folder_name[i]<-"interaction_center"
  }
}

v_teoretical<-c(126,151,269,144,322,272,20,23,27)

df_select<-df_select%>%mutate(only_docking=F)
df_select<-df_select%>%mutate(only_teoretical=F)
df_select<-df_select%>%mutate(both=F)

df_select<-df_select%>%mutate(only_docking_resno=F)
df_select<-df_select%>%mutate(only_teoretical_resno=F)
df_select<-df_select%>%mutate(both_resno=F)
j<-1
for (j in 1:nrow(df_select)) {
  df_seq<-read.csv(paste0(part_start,"MD_analysis/fin_data/str_data/",df_select$receptor[j],".csv"),stringsAsFactors = F)
  df_selected<-df_merge%>%filter(receptor==df_select$receptor[j] )
  df_selected<-df_selected%>%filter(ligand==df_select$ligand[j] )
  df_inteteraction<-read.csv(paste0('din/',df_selected$folder_name[1],"/", df_selected$name.x[1],".csv"),stringsAsFactors = F)
  df_inteteraction<-df_inteteraction%>%select(resno,resid,persent_interactions)

  df_inteteraction<-unique(df_inteteraction)
  if (nrow(df_selected)>1){
    for (i in 2:nrow(df_selected)) {
      df_inteteraction_add<-read.csv(paste0('din/',df_selected$folder_name[i],"/", df_selected$name.x[i],".csv"),stringsAsFactors = F)
      df_inteteraction_add<-df_inteteraction_add%>%select(resno,resid,persent_interactions)
      df_inteteraction_add<-unique(df_inteteraction_add)
      df_inteteraction<-rbind(df_inteteraction,df_inteteraction_add)
    }
  }
  df_inteteraction<-df_inteteraction%>%group_by(resno)%>%
    mutate(max_persent_interactions=max(persent_interactions))%>%
    filter(max_persent_interactions==persent_interactions)
  
  df_inteteraction<-unique(df_inteteraction)
  df_inteteraction<-df_inteteraction%>%select(resno,resid,persent_interactions)
  df_inteteraction<-df_inteteraction%>%mutate(experiment=F)
  df_inteteraction$experiment[df_inteteraction$resno%in%v_teoretical]<-T
  df_inteteraction<-df_inteteraction%>%mutate(docking=F)
  df_inteteraction$docking[df_inteteraction$persent_interactions==100]<-T

  df_seq<-left_join(df_seq,df_inteteraction,by=c("resno", "resid"))
  write.csv(df_seq,paste0(part_start,"MD_analysis/fin_data/docking_data/",df_select$receptor[j],"_",df_select$ligand[j],".csv"),row.names = F)
  df_seq<-df_seq%>%mutate(amino=paste0(resno,resid))
  df_docking_only<-df_seq[!df_seq$experiment&df_seq$docking,]
  df_experiment_only<-df_seq[df_seq$experiment&!df_seq$docking,]
  df_both<-df_seq[df_seq$experiment&df_seq$docking,]

  df_select$only_docking[j]<-paste0(df_docking_only$amino,collapse = " ")
  df_select$only_teoretical[j]<-paste0(df_experiment_only$amino,collapse = " ")
  df_select$both[j]<-paste0(df_both$amino,collapse = " ")
  df_select$only_docking_resno[j]<-paste0(df_docking_only$resno,collapse = " ")
  df_select$only_teoretical_resno[j]<-paste0(df_experiment_only$resno,collapse = " ")
  df_select$both_resno[j]<-paste0(df_both$resno,collapse = " ")
 
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('cd ', part_name,"\n\n",
                      'color Display Background white\n',
                      'color Labels Atoms black\n',
                      'color Labels Bonds black\n\n',
                      'mol new {',paste0("din/tost/",
                                         
                                         df_selected$receptor[1],"_",
                                         df_selected$ligand[1],"/",
                                         df_selected$structure_order[1],"_",
                                         df_selected$name.x[1]),'} type {pdb}')
  for (i in 2:nrow(df_selected)) {
    df_tcl[i,1]<-paste0('mol addfile  {',paste0("din/tost/",
                                                
                                                df_selected$receptor[i],"_",
                                                df_selected$ligand[i],"/",
                                                df_selected$structure_order[i],"_",
                                                df_selected$name.x[i]),'} type {pdb}')
    
  }
  df_tcl[i+1,1]<-paste0('mol modselect 0 0 chain A\n',
                        'mol modmaterial 0 0 Transparent\n',
                        'mol modstyle 0 0 NewCartoon\n',
                        'mol modselect 1 0 chain B\n',
                        'mol addrep 0\n',
                        'mol modselect 1 0 chain B\n',
                        'mol modstyle 1 0 Licorice\n',
                        'mol addrep 0\n',
                        'mol modselect 2 0 chain A and resid ',
                        df_select$only_docking_resno[j],'\n',
                        'mol modcolor 2 0 ColorID 0\n',
                        'mol modstyle 2 0 NewCartoon\n',
                        'mol addrep 0\n',
                        'mol modselect 3 0 chain A and resid ',
                        df_select$only_teoretical_resno[j],'\n',
                        'mol modmaterial 3 0 Transparent\n',
                        'mol modcolor 3 0 ColorID 1\n',
                        'mol modstyle 3 0 NewCartoon\n',
                        'mol addrep 0\n',
                        'mol modselect 4 0 chain A and resid ',
                        df_select$both_resno[j],'\n',
                        'mol modcolor 4 0 ColorID 1\n',
                        'mol modstyle 4 0 NewCartoon\n')#,
  
  write.table(df_tcl,paste0("din/tcl/",df_select$receptor[j],"_",df_select$ligand[j],".tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
}
write.csv(df_select,paste0("all_interactions_plot.csv"),row.names = F)
