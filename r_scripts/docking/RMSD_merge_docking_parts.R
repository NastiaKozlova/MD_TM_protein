part_analysis <- commandArgs(trailingOnly=TRUE)
#group ligand structures
#group ligand structures
library(bio3d)
library(dplyr)
library(ggplot2)
v_rmsd<-3.5

print(Sys.time())

setwd(part_analysis)
df_all<-read.csv(paste0(part_analysis,"df_separate.csv"),stringsAsFactors = F)
#df_all<-df_all%>%filter(surf=="surf")
#df_all<-unique(df_all)
df_analysis<-df_all%>%select(receptor,ligand)
df_analysis<-unique(df_analysis)
df_analysis<-df_analysis%>%mutate(receptor_ligand=paste0(receptor,"_",ligand))

part<-paste0(part_analysis,"din/")
setwd(part)
if(dir.exists(paste0(part,"fin_merged"))) {system(command = paste0("rm -r ",part,"fin_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"structure_merged"))) {system(command = paste0("rm -r ",part,"structure_merged"),ignore.stdout=T,wait = T)}
if(dir.exists(paste0(part,"groups_merged"))) {system(command = paste0("rm -r ",part,"groups_merged"),ignore.stdout=T,wait = T)}
if (!dir.exists("RMSD_merged")) {dir.create("RMSD_merged")}
if (!dir.exists("groups_merged")) {dir.create("groups_merged")}
if (!dir.exists("structure_merged")) {dir.create("structure_merged")}
if (!dir.exists("fin_merged")) {dir.create("fin_merged")}

df_structure<-read.csv("RMSD_group.csv",stringsAsFactors = F)

df_structure<-df_structure%>%select(models.x, ligand_center,receptor, ligand, center,
                                    surf, x, y,z)
df_structure<-unique(df_structure)
df_structure<-df_structure%>%filter(surf=="surf")
#df_structure<-left_join(df_structure,df_all,by = join_by(receptor, ligand, center))

q<-1
i<-1
print(Sys.time())

for (q in 1:nrow(df_analysis)) {
  if(!file.exists(paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"))){
    df_structure_TEMP<-df_structure%>%filter(receptor==df_analysis$receptor[q])
    df_structure_TEMP<-df_structure_TEMP%>%filter(ligand==df_analysis$ligand[q])
    
    df_structure_merge<-left_join(df_structure_TEMP,df_structure_TEMP,by=c("receptor", "ligand","surf"),
                                  relationship = "many-to-many")
    df_structure_merge<-df_structure_merge%>%mutate(RMSD=NA)

    df_structure_merge<-df_structure_merge%>%mutate(x=abs(x.x-x.y))
    df_structure_merge<-df_structure_merge%>%mutate(y=abs(y.x-y.y))
    df_structure_merge<-df_structure_merge%>%mutate(z=abs(z.x-z.y))
    print(nrow(df_structure_merge))
    df_structure_merge<-df_structure_merge%>%filter(x<=30)
    df_structure_merge<-df_structure_merge%>%filter(y<=30)
    df_structure_merge<-df_structure_merge%>%filter(z<=30)
    
    df_structure_merge<-df_structure_merge%>%filter(x.x<=x.y)
    df_structure_merge<-df_structure_merge%>%filter(y.x<=y.y)
    df_structure_merge<-df_structure_merge%>%filter(z.x<=z.y)
    
    print(nrow(df_structure_merge)) 
    print(paste(df_analysis$receptor_ligand[q],Sys.time()))
    
#    df_structure_merge<-df_structure_merge%>%mutate(name.x=paste0("str_fin/",df_structure_merge$receptor[j],"_",
#                                                                  df_structure_merge$ligand[j],"_",
#                                                                  df_structure_merge$center.x[j],"_",
#                                                                  df_structure_merge$models.x.x[j]))
#    df_structure_merge<-df_structure_merge%>%mutate(name.y=paste0("str_fin/",df_structure_merge$receptor[j],"_",
#                                                                  df_structure_merge$ligand[j],"_",
#                                                                  df_structure_merge$center.y[j],"_",
#                                                                  df_structure_merge$models.x.y[j]))
#    df_structure_merge<-df_structure_merge%>%mutate(test.x=file.exists(df_structure_merge$name.x))
#    df_structure_merge<-df_structure_merge%>%mutate(test.y=file.exists(df_structure_merge$name.y))
#    df_structure_merge_x<-df_structure_merge%>%filter(!test.x)
#    df_structure_merge_y<-df_structure_merge%>%filter(!test.y)
    if(nrow(df_structure_merge)>0){
      for (j in 1:nrow(df_structure_merge)) {
        pdb_1<-read.pdb(paste0("str_fin/",df_structure_merge$receptor[j],"_",
                               df_structure_merge$ligand[j],"_",
                               df_structure_merge$center.x[j],"_",
                               df_structure_merge$models.x.x[j]))
        pdb_2<-read.pdb(paste0("str_fin/",df_structure_merge$receptor[j],"_",
                               df_structure_merge$ligand[j],"_",
                               df_structure_merge$center.y[j],"_",
                               df_structure_merge$models.x.y[j]))
        
        df_structure_merge$RMSD[j]<-rmsd(pdb_1,pdb_2)
      }
#      print(paste(df_analysis$receptor_ligand[q],nrow(df_structure_merge_test)))
      
      df_structure_merge<-df_structure_merge%>%filter(RMSD<50)
      
      write.csv(df_structure_merge,paste0("RMSD_merged/",df_analysis$receptor_ligand[q],".csv"),row.names=F)
#      df_structure_merge_test<-df_structure_merge%>%filter(is.na(RMSD))

      print(paste(df_analysis$receptor_ligand[q],Sys.time()))
    }
  }
}
print(Sys.time())
