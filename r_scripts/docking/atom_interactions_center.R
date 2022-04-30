part_start <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
setwd(paste0(part_start,"MD_analysis/docking/docking_first/"))
df_all<-read.csv(paste0("df_merge_center.csv"),stringsAsFactors = F)
#df_all<-df_all%>%select(receptor,ligand,center)
#df_all<-unique(df_all)

i<-18
if(!dir.exists("complex_structure_center")){dir.create("complex_structure_center")}
if(!dir.exists("plot_complex_tcl")){dir.create("plot_complex_tcl")}
df_merge<-df_all
df_merge<-df_merge%>%mutate(complex_name=paste0(receptor,"_",ligand,"_",center,"_",grop_number,"_",models.x))
df_merge<-df_merge%>%mutate(pdb_name=paste0(receptor,"_",ligand,"_",center,"_",size_of_group))
df_merge$models.y<-NULL
df_merge$RMSD<-NULL
df_merge<-unique(df_merge)
i<-1
for (i in 1:nrow(df_merge)) {
  df_interactions<-read.csv(paste0("din/interaction_complex/",df_merge$receptor[i],"_",df_merge$ligand[i],"_",df_merge$center[i],".csv"),stringsAsFactors = F)
  df_interactions<-df_interactions%>%filter(total_persent_interactions>0)
  df_interactions<-df_interactions%>%select(resid,resno,receptor,ligand,size_of_group,total_persent_interactions)
  df_interactions<-unique(df_interactions)
  df_merge_TEMP<-df_merge[i,]
  df_interactions_TEMP<-semi_join(df_interactions,df_merge_TEMP,by = c("receptor", "ligand", "size_of_group"))
  
  receptor_name<-paste0("receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0("din/str_fin/",df_merge$complex_name[i])
  
  #  complex_name<-paste0("complex_structure_center/",df_merge$complex_name[i],".pdb")
  #  ligand_name<-paste0(part_name,"str_fin/",df_merge$name.x[i])
  receptor<-read.pdb(receptor_name)
  ligand<-read.pdb(ligand_name)
  df_protein<-protein$atom
  df_ligand<-ligand$atom
  protein_binding<-binding.site(protein,ligand,cutoff = 12)
  
  if(length(protein_binding)>0){
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
    df_interaction<-df_test%>%group_by(resno.x)%>%mutate(length_test=min(length))
    df_interaction<-df_interaction%>%group_by(resno.x)%>%filter(length_test==length)
    
    df_interaction<-ungroup(df_interaction)
    
    df_TEMP<-df_interaction[df_interaction$resno.x%in%df_interactions_TEMP$resno,]
    
    if (nrow(df_TEMP)==0){
      print(paste0(i,nrow(df_TEMP),nrow(df_interaction),collapse = " "))
    }
  }
} 
print(nrow(df_merge))
df_merge<-df_merge%>%filter(!is.na(models.x))
print(nrow(df_merge))
i<-1
setwd("din")
for (i in 1:nrow(df_merge)) {
  df_hbonds<-read.csv(paste0(part_start,"MD_analysis/din/",df_merge$receptor[i],"/hbonds_8.csv"),stringsAsFactors = F)
  df_hbonds<-df_hbonds%>%filter(persent>50)
  df_interactions<-read.csv(paste0("interaction_complex/",df_merge$receptor[i],"_",df_merge$ligand[i],"_",df_merge$center[i],".csv"),stringsAsFactors = F)
  df_interactions<-df_interactions%>%filter(total_persent_interactions>0)
  df_interactions<-df_interactions%>%select(resid,resno,receptor,ligand,size_of_group,total_persent_interactions)
  df_interactions<-unique(df_interactions)
  df_merge_TEMP<-df_merge[i,]
  df_interactions_TEMP<-semi_join(df_interactions,df_merge_TEMP,by = c("receptor", "ligand", "size_of_group"))
  receptor_name<-paste0(part_start,"MD_analysis/docking/docking_first/receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0("str_fin/",df_merge$complex_name[i])
  protein<-read.pdb(receptor_name)
  ligand<-read.pdb(ligand_name)
  v_binding_test<-binding.site(protein,ligand,cutoff = 12)$resno
  df_protein<-protein$atom
  df_ligand<-ligand$atom
  df_protein<-df_protein[df_protein$resno%in%df_interactions_TEMP$resno,]
  df_protein<-unique(df_protein[df_protein$resno%in%v_binding_test,])
  if(nrow(df_protein)>0){
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
    df_test<-df_test%>%filter(alt.x!=alt.y)
    df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                              eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
    df_interaction<-df_test%>%group_by(resno.x)%>%mutate(length_test=min(length))
    df_interaction<-df_interaction%>%group_by(resno.x)%>%filter(length_test==length)
    df_interaction<-ungroup(df_interaction)
    
    df_interaction<-df_interaction%>%group_by(eleno.y)%>%mutate(length_test=min(length))
    df_interaction<-df_interaction%>%group_by(eleno.y)%>%filter(length_test==length)
    df_interaction<-ungroup(df_interaction)
    
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0('cd ', part_name,"din/complex_structure_center/\n\n",
                        'mol new {',df_merge$pdb_name[i],'.pdb} type {pdb}')
    b<-paste('(resid ',df_interaction$resno.x,' and name ',df_interaction$elety.x," and resname ",df_interaction$resid.x,")")
    a<-paste0(b,collapse = " or ")
    df_tcl[1,2]<-paste0('set all [atomselect ',(i-1),' "',a,'"]\n',
                        'set i ',(i-1),"\n",
                        'foreach atom [$all list] {\n',
                        '  label add Atoms ',(i-1),'/$atom\n',
                        '  incr i\n}\n',
                        '$all delete\n\n',
                        'color Display Background white\n',
                        'color Labels Atoms black\n',
                        'color Labels Bonds black\n\n')
    
    df_tcl[1,3]<-paste0('mol modselect 0 ',i-1,' protein\n',
                        'mol modstyle 0 ' ,i-1, ' NewCartoon\n',
                        'mol selection resid ',paste0(unique(df_hbonds$number),collapse = " "))
    df_tcl[1,4]<-paste0('mol material Opaque\n',
                        'mol addrep ',(i-1),'\n',
                        'mol modstyle 1 ',(i-1),' Licorice\n',
                        'mol modcolor 1 ',(i-1),' ColorID 0\n')
    if (nrow(df_interaction)>0){
      df_tcl[1,5]<-paste0('mol selection resname ',paste0(unique(df_interaction$resid.y),collapse = " "))
      df_tcl[1,6]<-paste0('mol material Opaque\n',
                          'mol addrep ',(i-1),'\n',
                          'mol modstyle 2 ',(i-1),' Licorice\n',
                          'mol modcolor 2 ',(i-1),' ColorID 16\n')
      df_tcl[1,7]<-paste0('mol selection (resid ',paste0(unique(df_interaction$resno.x),collapse = " "),")")
      df_tcl[1,8]<-paste0(' mol material Opaque\n',
                          'mol addrep ',(i-1),'\n',
                          'mol modstyle 3 ',(i-1),' Licorice\n',
                          'mol modcolor 3 ',(i-1),' ColorID 1')
      
      for (p in 1:nrow(df_interaction)) {
        df_tcl[(p+1),1]<-paste0('set atomID1 [[atomselect ',(i-1),' "(resid ',df_interaction$resno.x[p],
                                ' and name ',df_interaction$elety.x[p],
                                " and resname ",df_interaction$resid.x[p],')"] list]')
        df_tcl[(p+1),2]<-paste0('set atomID2 [[atomselect ',(i-1),
                                ' "(x > ',df_interaction$x.y[p]-0.5,' and x < ',df_interaction$x.y[p]+0.5,
                                ' and y > ',df_interaction$y.y[p]-0.5,' and y < ',df_interaction$y.y[p]+0.5,
                                ' and z > ',df_interaction$z.y[p]-0.5,' and z < ',df_interaction$z.y[p]+0.5,')"] list]')
        
        df_tcl[(p+1),3]<-paste0('label add Bonds ',(i-1),'/$atomID1 ',(i-1),'/$atomID2')
      }
      df_tcl[is.na(df_tcl)]<-""
      write.csv(df_tcl,paste0("make_picture_center_tcl/",df_merge$complex_name[i],"_test.tcl"),row.names = F)
    }
  }
} 

df_tcl<-read.csv(paste0("make_picture_center_tcl/",df_merge$complex_name[1],"_test.tcl"),stringsAsFactors = F)
i<-2
for (i in 2:nrow(df_merge)) {
  df_tcl_add<-read.csv(paste0("make_picture_center_tcl/",df_merge$complex_name[i],"_test.tcl"),stringsAsFactors = F)
  df_tcl<-rbind(df_tcl,df_tcl_add)
}
write.table(df_tcl,paste0("make_picture_center_tcl.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
