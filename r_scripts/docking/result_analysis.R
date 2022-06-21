
part_analysis <- commandArgs(trailingOnly=TRUE)
part_analysis<-"/media/nastia/SSD_Alya/NaPi2b/MD_TM_protein/MD_analysis/docking/docking_first/"
library(ggplot2)
library(dplyr)
library(rstatix)
a<-seq(from=0,to=10000,by=100)
setwd(part_analysis)
setwd("din")
df_structure_RMSD<-read.csv("df_merge_structure_log.csv",stringsAsFactors =  F)

#df_structure_RMSD
df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(group_size=n())
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-df_structure_RMSD%>%select(name.x,group_size,receptor,ligand)
df_structure_RMSD<-unique(df_structure_RMSD)
print(nrow(df_structure_RMSD))

df_structure_RMSD<-df_structure_RMSD%>%group_by(ligand)%>%mutate(median_group_size=median(group_size))
df_structure_median<-df_structure_RMSD%>%select(receptor,ligand,median_group_size)
df_structure_median<-ungroup(df_structure_median)
df_structure_median<-unique(df_structure_median)
p<-ggplot(data=df_structure_RMSD)+
  geom_freqpoly(aes(x=group_size),binwidth=10)+
  
  geom_segment(aes(x = median_group_size, y = -Inf, xend = median_group_size, yend = Inf),data=df_structure_median)+
  facet_grid(ligand~receptor)+
  scale_x_continuous(breaks=a,labels=a)+
  theme_bw()
p
df_structure<-df_structure_RMSD%>%filter(group_size>median_group_size)
df_structure_RMSD<-read.csv("df_merge_structure_log.csv",stringsAsFactors =  F)
df_structure_RMSD<-df_structure_RMSD[df_structure_RMSD$name.x%in%df_structure$name.x,]
df_structure_RMSD<-df_structure_RMSD%>%group_by(ligand)%>%mutate(ligand_size=n())
data_down_results_pwc <- df_structure_RMSD %>%
  kruskal_test(affinity ~ ligand)
#df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(group_size=n())
#df_structure_RMSD<-df_structure_RMSD%>%group_by(name.x)%>%mutate(group_size>400)
data_down_results_pwc <- df_structure_RMSD %>%
  pairwise_wilcox_test(affinity ~ ligand, p.adjust.method = "BH")
#df_structure_RMSD<-df_structure_RMSD%>%mutate(name=paste(receptor,ligand,center,sep = "_"))
#df_tost<-df_structure_RMSD %>%
#  group_by(name) %>%
#  wilcox_test(affinity ~ name.x)


df_structure_RMSD<-unique(df_structure_RMSD)
print(nrow(df_structure_RMSD))
b<-seq(from=-5,to=5,by=0.1)
p<-ggplot(data=df_structure_RMSD)+
  geom_boxplot(aes(y=affinity,x=ligand,colour=as.factor(group_size)))+
  scale_y_continuous(breaks = b,labels = b)

  #  geom_segment(aes(x = mediana, y = -Inf, xend = mediana, yend = Inf),data=df_structure_median)+
    facet_grid(ligand~receptor)+
  #  scale_x_continuous(breaks=z,labels=z)+
  theme_bw()
p
df_structure_RMSD<-df_structure_RMSD%>%filter(group_size>100)
print(nrow(df_structure_RMSD))


df_structure_RMSD<-df_structure_RMSD%>%select(group_size)

df_structure_RMSD<-df_structure_RMSD%>%group_by(group_size)%>%mutate(number=n())

df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_RMSD<-unique(df_structure_RMSD)

p<-ggplot(data=df_structure_RMSD)+
  geom_line(aes(x=group_size,y=number))+
  geom_point(aes(x=group_size,y=number))+
#  geom_segment(aes(x = mediana, y = -Inf, xend = mediana, yend = Inf),data=df_structure_median)+
#  facet_grid(ligand~receptor)+
#  scale_x_continuous(breaks=z,labels=z)+
  theme_bw()

df_structure_RMSD$name.x<-NULL
df_structure_RMSD<-unique(df_structure_RMSD)

#df_structure_RMSD_add<-data.frame(matrix(ncol=ncol(df_structure_RMSD),nrow = max(df_structure_RMSD$group_size)))
#colnames(df_structure_RMSD_add)<-colnames(df_structure_RMSD)
#df_structure_RMSD_add$number<-0
#df_structure_RMSD_add$group_size<-1:nrow(df_structure_RMSD_add)
#print(nrow(df_structure_RMSD_add))
#df_structure_RMSD_add<-df_structure_RMSD_add[!df_structure_RMSD_add$group_size%in%df_structure_RMSD$group_size,]
#print(nrow(df_structure_RMSD_add))
#df_structure_RMSD<-rbind(df_structure_RMSD,df_structure_RMSD_add)
df_structure_RMSD<-df_structure_RMSD%>%arrange(group_size)
p<-30
df_structure_RMSD<-df_structure_RMSD%>%mutate(mean_number=NA)
df_structure_RMSD<-df_structure_RMSD%>%mutate(mean_group_size=NA)
for (j in (1+p):(nrow(df_structure_RMSD)-p)) {
  df_structure_RMSD$mean_group_size[j]<-mean(df_structure_RMSD$group_size[(j-p):(j+p)],na.rm = F)
  df_structure_RMSD$mean_number[j]<-mean(df_structure_RMSD$number[(j-p):(j+p)],na.rm = F)
}
df_structure_RMSD<-df_structure_RMSD%>%mutate(ligand_receptor=paste(ligand,receptor))
df_structure_RMSD<-df_structure_RMSD%>%group_by(ligand_receptor)%>%mutate(mediana=median(group_size))
df_structure_RMSD<-ungroup(df_structure_RMSD)
df_structure_median<-df_structure_RMSD%>%select(receptor,ligand,ligand_receptor,mediana)
df_structure_median<-unique(df_structure_median)
z<-seq(from=0,to=1000,by=50)
p<-ggplot(data=df_structure_RMSD)+geom_point(aes(x=group_size,y=number))+
  geom_segment(aes(x = mediana, y = -Inf, xend = mediana, yend = Inf),data=df_structure_median)+
  facet_grid(ligand~receptor)+
  scale_x_continuous(breaks=z,labels=z)+
  theme_bw()
#+guides(color = "none", size = "none")
ggsave(p,filename = paste0("group_size_ligand_receptor_center.png"), width = 20, height = 20, units = c("cm"), dpi = 200 )
data_down_results_pwc <- df_structure_RMSD %>%
  pairwise_wilcox_test(affinity ~ name, p.adjust.method = "BH")