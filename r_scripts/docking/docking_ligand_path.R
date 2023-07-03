part_start <- commandArgs(trailingOnly=TRUE)
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/")
library(dplyr)
library(ggplot2)
library(bio3d)
setwd(part_name)
df_merge<-read.csv(file = paste0(part_name,"din/df_structure.csv"),stringsAsFactors = F)



if (!dir.exists(paste0(part_name,"tost/"))) {dir.create(paste0(part_name,"tost/"))}
if (!dir.exists(paste0(part_name,"plot/"))) {dir.create(paste0(part_name,"plot/"))}
if (!dir.exists(paste0(part_name,"all_interactions/"))) {dir.create(paste0(part_name,"all_interactions/"))}
#"all_interactions/"
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
i<-1
j<-1
df_inteteraction<-read.csv(paste0('din/',df_selected$folder_name[1],"/", df_merge$name.x[1],".csv"),stringsAsFactors = F)

df_inteteraction<-df_inteteraction%>%filter(persent_interactions==100)
df_inteteraction<-df_inteteraction%>%mutate(receptor=df_merge$receptor[1])
df_inteteraction<-df_inteteraction%>%mutate(ligand=df_merge$ligand[1])
df_inteteraction<-df_inteteraction%>%mutate(name.x=df_merge$name.x[1])
if (nrow(df_merge)>1){
  for (i in 2:nrow(df_merge)) {
  df_inteteraction_add<-read.csv(paste0('din/',df_selected$folder_name[1],"/", df_merge$name.x[1],".csv"),stringsAsFactors = F)
  
  df_inteteraction_add<-df_inteteraction_add%>%filter(persent_interactions==100)
  df_inteteraction_add<-df_inteteraction_add%>%mutate(receptor=df_merge$receptor[i])
  df_inteteraction_add<-df_inteteraction_add%>%mutate(ligand=df_merge$ligand[i])
  df_inteteraction_add<-df_inteteraction_add%>%mutate(name.x=df_merge$name.x[i])
  df_inteteraction<-rbind(df_inteteraction,df_inteteraction_add)
  }
}
i<-1
for (i in 1:nrow(df_select)) {
  df_temp<-df_merge%>%filter(receptor==df_select$receptor[i])
  df_temp<-df_temp%>%filter(ligand==df_select$ligand[i])
#  df_inteteraction<-read.csv(paste0('din/',df_selected$folder_name[1],"/", df_selected$name.x[1],".csv"),stringsAsFactors = F)
  df_interacting<-df_inteteraction[df_inteteraction$name.x%in%c(df_temp$name.x),]
  df_interacting<-df_interacting%>%select(resno, resid,receptor, ligand)
  df_interacting<-unique(df_interacting)
 
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('cd ', part_name,"\n\n",
                      'color Display Background white\n',
                      'color Labels Atoms black\n',
                      'color Labels Bonds black\n\n',
                      'mol new {',paste0("din/tost/",
                                         
                                         df_temp$receptor[1],"_",
                                         df_temp$ligand[1],"/",
                                         df_temp$structure_order[1],"_",
                                         df_temp$name.x[1]),'} type {pdb}')
  for (j in 2:nrow(df_temp)) {
    df_tcl[j,1]<-paste0('mol addfile  {',paste0("din/tost/",
                                                
                                                df_temp$receptor[j],"_",
                                                df_temp$ligand[j],"/",
                                                df_temp$structure_order[j],"_",
                                                df_temp$name.x[j]),'} type {pdb}')
    
  }
  df_tcl[j+1,1]<-paste0('mol modselect 0 ',i-1,' chain A\n',
                        'mol modmaterial 0 ',(i-1),' Transparent\n',
                        'mol modstyle 0 ' ,i-1, ' NewCartoon\n',
                        'mol modselect 1 ',i-1,' chain B\n',
                        'mol addrep ',(i-1),'\n',
                        'mol modselect 1 ',i-1,' chain B\n',
                        'mol modstyle 1 ' ,i-1, ' Licorice\n',
                        'mol addrep ',(i-1),'\n',
                        'mol modselect 2 ',i-1,' chain A and resid ',
                        
                        paste0(df_interacting$resno,collapse = " "),'\n',
                        #'mol addrep ',(i-1),'\n',
                        'mol modcolor 2 0 ColorID 1\n',
                        'mol modstyle 2 ' ,i-1, ' NewCartoon\n')#,
  
  write.table(df_tcl,paste0("din/tcl/",df_select$receptor[i],"_",df_select$ligand[i],".tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
#vmd
  write.csv(df_interacting,paste0("all_interactions/",df_select$receptor[i],"_",df_select$ligand[i],".csv"),row.names = F)
}
