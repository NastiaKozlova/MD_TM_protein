part_analysis <- commandArgs(trailingOnly=TRUE)

library(bio3d)
library(dplyr)
library(ggplot2)
part_hbonds<-strsplit(part_analysis,split = "/",fixed=T)[[1]]
part_topology<-paste(part_hbonds[1:(length(part_hbonds)-3)],collapse = "/")
part_hbonds<-paste(part_hbonds[1:(length(part_hbonds)-2)],collapse = "/")

part_hbonds<-paste0(part_hbonds,"/")
part_topology<-paste0(part_topology,"/start/")
df_topology<-read.csv(paste0(part_topology,"df_topology.csv"),stringsAsFactors = F)

df_all<-read.csv(paste0(part_analysis,"df_all.csv"),stringsAsFactors = F)
df_all<-df_all%>%select(receptor,ligand)
df_all<-unique(df_all)
part_name<-paste0(part_analysis,"din/")
setwd(part_name)
df_merge<-read.csv(paste0(part_name,"df_merge_structure_log_center.csv"),stringsAsFactors = F)
df_merge<-df_merge%>%select(name.x,receptor,ligand, size_of_group)
df_merge<-unique(df_merge)

if (dir.exists(paste0("make_picture_tcl_center/"))) {system(command = paste0("rm -r ",part_analysis,"din/make_picture_tcl_center/"))}
if(!dir.exists("make_picture_tcl_center")){dir.create("make_picture_tcl_center")}
df_merge<-df_merge%>%mutate(complex_name=paste0(receptor,"_",ligand,"_",size_of_group))
i<-1
j<-1
df_merge<-df_merge%>%mutate(exists=F)
df_merge<-df_merge%>%filter(receptor=="charmm-gui-1717818438")
for (i in 1:nrow(df_merge)) {
  if(file.exists(paste0("interaction_center/",df_merge$name.x[i],".csv"))){
    if(file.exists(paste0("complex_structure_center/",df_merge$name.x[i]))){
      df_merge$exists<-T
    }
  }
}
i<-1
df_merge<-df_merge%>%filter(exists)
for (i in 1:nrow(df_merge)) {
  df_topology_seclected<-df_topology%>%mutate(selected=F)
  
  df_hbonds<-read.csv(paste0(part_hbonds,"fin_data/str_data/",df_merge$receptor[i],".csv"),stringsAsFactors = F)
  df_hbonds<-df_hbonds%>%filter(hbonds>50)
  df_hbonds<-df_hbonds%>%filter(abs(z)<20)
  df_interactions<-read.csv(paste0("interaction_center/",df_merge$name.x[i],".csv"),stringsAsFactors = F)
  df_interactions<-df_interactions%>%filter(persent_interactions==100)
  df_interactions<-df_interactions%>%select(resid,resno,persent_interactions)
  df_interactions<-unique(df_interactions)
  
  for (p in 1:nrow(df_topology_seclected)) {
    v_domain<-df_interactions$resno[df_interactions$resno%in%c(df_topology_seclected$seq_beg[p]:df_topology_seclected$seq_end[p])]
    if(length(v_domain)>0){
      df_topology_seclected$selected[p]<-T
    }
  }
  df_topology_seclected<-df_topology_seclected%>%filter(selected)
  
  receptor_name<-paste0(part_analysis,"receptor_start/",df_merge$receptor[i],".pdb")
  ligand_name<-paste0(part_name,"str_fin/",df_merge$name.x[i])
  protein<-read.pdb(receptor_name)
  ligand<-read.pdb(ligand_name)
  v_binding_test<-binding.site(protein,ligand,cutoff = 12)$resno
  df_protein<-protein$atom
  df_ligand<-ligand$atom
  df_protein<-df_protein[df_protein$resno%in%df_interactions$resno,]
  df_protein<-unique(df_protein[df_protein$resno%in%v_binding_test,])
  
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
  df_test<-df_test%>%filter(length<5)
  df_test<-df_test%>%filter(alt.x!=alt.y)
  df_test<-df_test%>%select(eleno.x,  elety.x,  alt.x,    resid.x,  resno.x,x.x,y.x,z.x,
                            eleno.y,  elety.y,  alt.y,    resid.y,  resno.y,x.y,y.y,z.y, length)
  df_interaction<-df_test%>%group_by(resno.x)%>%mutate(length_test=min(length))
  df_interaction<-df_interaction%>%group_by(resno.x)%>%filter(length_test==length)
  df_interaction<-ungroup(df_interaction)
  v_selected<-c()
  for (q in 1:nrow(df_topology_seclected)) {
    v_selected<-c(v_selected,df_topology_seclected$seq_beg[q]:df_topology_seclected$seq_end[q])
  }
  v_selected<-paste(v_selected,collapse = " ")
  
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
  df_tcl[1,1]<-paste0('cd ', part_name,"complex_structure_center/\n\n",
                      'mol new {',df_merge$name.x[i],'} type {pdb}')
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
  
  df_tcl[1,3]<-#paste0('mol modselect 0 ',i-1,' protein\n',
               paste0('mol modselect 0 ',i-1,' protein and resid ',v_selected,'\n',
                      'mol modmaterial 0 ',(i-1),' Transparent\n',
                      'mol modstyle 0 ' ,i-1, ' NewCartoon\n')#,
  if (nrow(df_interaction)>0){
    df_tcl[1,4]<-paste0('mol selection resname ',paste0(unique(df_interaction$resid.y),collapse = " "))
    df_tcl[1,5]<-paste0('mol modmaterial 1 ',(i-1),' Opaque\n',
                        'mol addrep ',(i-1),'\n',
                        'mol modstyle 1 ',(i-1),' Licorice \n',
                        'mol modcolor 1 ',(i-1),' Name\n')
    df_tcl[1,6]<-paste0('mol selection (resid ',paste0(unique(df_interaction$resno.x),collapse = " "),")")
    df_tcl[1,7]<-paste0(' mol modmaterial 2 ',(i-1),' Opaque\n',
                        'mol addrep ',(i-1),'\n',
                        'mol modstyle 2 ',(i-1),' Licorice 0.1 12 12 \n',
                        'mol modcolor 2 ',(i-1),' Type')
    
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
  }
  v_selected<-c()
  for (q in 1:nrow(df_topology_seclected)) {
    v_selected<-c(v_selected,df_topology_seclected$seq_beg[q]:df_topology_seclected$seq_end[q])
  }
  v_selected<-paste(v_selected,collapse = " ")
 # df_tcl[p+2,1]<-paste0('mol modselect 0 ',i-1,' protein and resid ',v_selected,'\n',
#                        'mol modmaterial 0 ',(i-1),' Opaque\n',
#                        'mol modstyle 0 ' ,i-1, ' NewCartoon\n')#,
#  df_TMD25<-paste(df_topology[df_topology$type%in%c("TM5","TM2"),]
  v_loop1<-paste(df_topology$seq_beg[df_topology$type%in%c("loop1")]:df_topology$seq_end[df_topology$type%in%c("loop1")],collapse = " ")
  v_loop2<-paste(df_topology$seq_beg[df_topology$type%in%c("loop2")]:df_topology$seq_end[df_topology$type%in%c("loop2")],collapse = " ")
  v_TMD2<-paste(df_topology$seq_beg[df_topology$type%in%c("TM2")]:df_topology$seq_end[df_topology$type%in%c("TM2")],collapse = " ")
  v_TMD5<-paste(df_topology$seq_beg[df_topology$type%in%c("TM5")]:df_topology$seq_end[df_topology$type%in%c("TM5")],collapse = " ")
#  v_TMD25<-paste(c(df_TMD25$seq_beg[1]:df_TMD25$seq_end[1],df_TMD25$seq_beg[2]:df_TMD25$seq_end[2]),collapse = " ")
  v_qsss<-paste(c(158:161,431:434),collapse = " ")
  
#  df_tcl[p+2,2]<-paste0('mol modselect 0 ',i-1,' protein and resid ',v_selected,'\n',
#                        'mol modmaterial 0 ',(i-1),' Transparent\n',
#                      #  'mol addrep ',(i-1),'\n',
#                        'mol modstyle 0 ',(i-1),' NewCartoon \n',
#                        'mol modcolor 0 ',(i-1),' Name\n')
  df_tcl[p+2,1]<-paste0('mol addrep ',(i-1),'\n',
                        'mol modselect 3 ',i-1,' protein and resid ',v_loop1,'\n',
                        'mol modmaterial 3 ',(i-1),' Opaque\n',
                        'mol modstyle 3 ',(i-1),' NewCartoon \n',
                        'mol modcolor 3 ',(i-1),' ColorID 4\n')#,
  
  df_tcl[p+2,2]<-paste0('mol addrep ',(i-1),'\n',
                        'mol modselect 4 ',i-1,' protein and resid ',v_loop2,'\n',
                        'mol modmaterial 4 ',(i-1),' Opaque\n',
                        'mol modstyle 4 ',(i-1),' NewCartoon \n',
                        'mol modcolor 4 ',(i-1),' ColorID 21\n')#,
  df_tcl[p+2,3]<-paste0('mol addrep ',(i-1),'\n',
                        'mol modselect 5 ',i-1,' protein and resid ',v_TMD2,'\n',
                        'mol modmaterial 5 ',(i-1),' Opaque\n',
                        'mol modstyle 5 ',(i-1),' NewCartoon \n',
                        'mol modcolor 5 ',(i-1),' ColorID 5\n')#,

  df_tcl[p+2,4]<-paste0('mol addrep ',(i-1),'\n',
                        'mol modselect 6 ',i-1,' protein and resid ',v_TMD5,'\n',
                        'mol modmaterial 6 ',(i-1),' Opaque\n',
                        'mol modstyle 6 ',(i-1),' NewCartoon \n',
                        'mol modcolor 6 ',(i-1),' ColorID 0\n')#,
  
  df_tcl[p+2,5]<-paste0('mol addrep ',(i-1),'\n',
                        'mol modselect 7 ',i-1,' protein and resid ',v_qsss,'\n',
                        'mol modmaterial 7 ',(i-1),' Opaque\n',
                        'mol modstyle 7 ',(i-1),' NewCartoon \n',
                        'mol modcolor 7 ',(i-1),' ColorID 16\n')#,
  #df_tcl[p+2,6]<-paste0('mol addrep ',(i-1),'\n',
  #                      'mol modselect 8 ',i-1,' protein and resid ',paste0(unique(df_hbonds$number),collapse = ' '),'\n',
  #                      'mol modmaterial 8 ',(i-1),' Transparent\n',
  #                      'mol modstyle 8 ',(i-1),' QuickSurf \n',
  #                      'mol modcolor 8 ',(i-1),' Name\n')#,
  

  df_tcl[is.na(df_tcl)]<-""
  write.csv(df_tcl,paste0("make_picture_tcl_center/",df_merge$name.x[i],".tcl"),row.names = F)
}
df_tcl<-read.csv(paste0("make_picture_tcl_center/",df_merge$name.x[1],".tcl"),stringsAsFactors = F)
i<-2
for (i in 2:nrow(df_merge)) {
  df_tcl_add<-read.csv(paste0("make_picture_tcl_center/",df_merge$name.x[i],".tcl"),stringsAsFactors = F)
  df_tcl<-rbind(df_tcl,df_tcl_add)
}
write.table(df_tcl,paste0("make_picture_tcl_center.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na="")
