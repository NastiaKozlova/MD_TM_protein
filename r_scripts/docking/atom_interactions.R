part_start <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
paste0(part_start)
df_all<-read.csv(paste0(part_start,"MD_analysis/docking/docking_first/df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%select(receptor,ligand)
df_all<-unique(df_all)
part_name<-paste0(part_start,"MD_analysis/docking/docking_first/din/")
setwd(part_name)
i<-1
for (i in 1:nrow(df_all)) {
  if(!file.exists(paste0("interaction_fin/",df_all$receptor[i],"_",df_all$ligand[i],".csv"))){
    df_all$receptor[i]<-NA
  }
}
df_all<-df_all%>%filter(!is.na(receptor))
df_merge<-read.csv(paste0(part_name,"df_merge_structure_log.csv"),stringsAsFactors = F)
df_merge<-semi_join(df_merge,df_all)
df_merge<-df_merge%>%select(name.x,receptor,ligand, size_of_group)
df_merge<-unique(df_merge)
i<-18
if(!dir.exists("complex_structure")){dir.create("complex_structure")}
if(!dir.exists("plot_tcl")){dir.create("plot_tcl")}
df_merge<-df_merge%>%mutate(complex_name=paste0(receptor,"_",ligand,"_",size_of_group))
for (i in 1:nrow(df_merge)) {
  
#  df_interactions<-read.csv(paste0("interaction_fin/",df_merge$receptor[i],"_",df_merge$ligand[i],".csv"),stringsAsFactors = F)
#  df_interactions<-df_interactions%>%filter(total_persent_interactions>0)
#  df_interactions<-df_interactions%>%select(resid,resno,receptor,ligand,size_of_group,total_persent_interactions)
#  df_interactions<-unique(df_interactions)
#  df_merge_TEMP<-df_merge[i,]
#  df_interactions_TEMP<-semi_join(df_interactions,df_merge_TEMP,by = c("receptor", "ligand", "size_of_group"))
  receptor_name<-paste0(part_start,"MD_analysis/docking/docking_first/receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_name,"str_fin/",df_merge$name.x[i])
  protein<-read.pdb(receptor_name)
  ligand<-read.pdb(ligand_name)
  df_protein<-protein$atom
  df_ligand<-ligand$atom
  df_protein<-df_protein[df_protein$resno%in%df_interactions_TEMP$resno,]
  
  q<-1
  for (q in 1:nrow(df_protein)) {
    df_protein$alt[q]<-strsplit(df_protein$elety[q],split = "",fixed = T)[[1]][1]
  }
  for (q in 1:nrow(df_ligand)) {
    df_ligand$alt[q]<-strsplit(df_ligand$elety[q],split = "",fixed = T)[[1]][1]
  }
  df_ligand<-df_ligand%>%filter(alt!="C")
  df_protein<-df_protein%>%filter(alt!="C")
  
  df_test<-full_join(df_protein,df_ligand,by="type")
  df_test<-df_test%>%mutate(length=sqrt((x.x-x.y)^2+(y.x-y.y)^2+(z.x-z.y)^2))
  df_test<-df_test%>%filter(length<12)
  if (nrow(df_test)>0){
    df_test<-df_test%>%filter(alt.x!=alt.y)
    df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                              eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
    df_interaction<-df_test%>%group_by(resno.x)%>%mutate(length_test=min(length))
    df_interaction<-df_interaction%>%group_by(resno.x)%>%filter(length_test==length)
    
    df_interaction<-ungroup(df_interaction)
    
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0('cd ', part_name,"complex_structure/\n\n")
    df_tcl[1,2]<-paste0('mol new {',df_merge$complex_name[i],'.pdb} type {pdb}')
    b<-paste('(resid ',df_interaction$resno.x,' and name ',df_interaction$elety.x," and resname ",df_interaction$resid.x,")")
    a<-paste0(b,collapse = " or ")
    df_tcl[1,3]<-paste0('set all [atomselect ',(i-1),' "',a,'"]')
    df_tcl[1,4]<-paste0('set i ',(i-1))
    df_tcl[1,5]<-paste0('foreach atom [$all list] {')
    df_tcl[1,6]<-paste0('  label add Atoms ',(i-1),'/$atom')
    df_tcl[1,7]<-paste0('  incr i\n}')
    df_tcl[1,8]<-paste0('$all delete')
    for (p in 1:nrow(df_interaction)) {
      
      
      df_tcl[(p+1),1]<-paste0('set atomID1 [[atomselect ',(i-1),' "(resid ',df_interaction$resno.x[p],
                              ' and name ',df_interaction$elety.x[p],
                              " and resname ",df_interaction$resid.x[p],')"] list]')
      df_tcl[(p+1),2]<-paste0('set atomID2 [[atomselect ',(i-1),
                              ' "(resid ',df_interaction$resno.y[p],
                              ' and name ',df_interaction$elety.y[p],
                              " and resname ",df_interaction$resid.y[p],
                              ' and x > ',df_interaction$x.y[p]-0.25,' and x < ',df_interaction$x.y[p]+0.25,
                              ' and y > ',df_interaction$y.y[p]-0.25,' and y < ',df_interaction$y.y[p]+0.25,
                              ' and z > ',df_interaction$z.y[p]-0.25,' and z < ',df_interaction$z.y[p]+0.25,')"] list]')
      
      df_tcl[(p+1),3]<-paste0('label add Bonds ',(i-1),'/$atomID1 ',(i-1),'/$atomID2')
    }
    df_tcl[is.na(df_tcl)]<-""
    write.table(df_tcl,paste0("make_picture_tcl/",df_merge$complex_name[i],"_test.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
#    print(paste0(df_merge$complex_name[i],' resid ',paste0(unique(df_interaction$resno.x),collapse = " ")))
  }else{
    print(paste(nrow(df_test),i,df_merge_TEMP$complex_name[1], file.exists(paste0(part_name,"complex_structure/",df_merge_TEMP$complex_name[1],".pdb"))))
  }
} 
