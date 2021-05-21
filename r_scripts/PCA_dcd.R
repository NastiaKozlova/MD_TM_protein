library(bio3d)
library(dplyr)
library(ggplot2)

setwd("/home/nastia/mem/NaPi2b/MD/SLC34A2/main_NaPi2b/NaPi2b_POPC/")
dcd <- read.dcd("namd/step7_production.dcd")
pdb <- read.pdb("namd/step5_charmm2namd.pdb")

ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd,
               fixed.inds=ca.inds$xyz,
               mobile.inds=ca.inds$xyz)

dim(xyz) == dim(dcd)

rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz])
plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)

hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD")
lines(density(rd), col="gray", lwd=3)

rf <- rmsf(xyz[,ca.inds$xyz])
plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")

pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=bwr.colors(nrow(xyz)) )

hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)
plot(pc, col=grps)

df_pc<-data.frame(pc$au)
df_pc<-df_pc%>%mutate(number=(1:nrow(df_pc)))

a<-seq(from=0,to=690,by=20)
p_PCA<-ggplot(data = df_pc)+
  labs(title=paste("PC1 (A)"), y="PC1 (A)", x="Residue Position") +
  geom_line(aes(x =number,y=X1,colour="PCA 1"))+
  geom_line(aes(x =number,y=X2,colour="PCA 2"))+
  geom_line(aes(x =number,y=X3,colour="PCA 3"))+
  
  scale_x_continuous(breaks = a,labels = a)+
  theme_bw()#+
  #geom_text(x=median(df_Energy$VdW), y=100,label=round(median(df_Energy$VdW),digits = 1))+
  #geom_vline(xintercept = median(df_Energy$VdW))
ggsave(p_PCA,filename = paste0("PCA_NaPi2b_POPC.png"), width = 40, height = 20, units = c("cm"), dpi = 200 ) 
write.csv(df_pc,"df_pc.csv",row.names = F)
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")

p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file="pc1.pdb")
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file="pc2.pdb")
p3 <- mktrj.pca(pc, pc=2,b=pc$au[,3], file="pc3.pdb")

write.ncdf(p1, "trj_pc1.nc")

cij<-dccm(xyz[,ca.inds$xyz])
plot(cij)

pymol.dccm(cij, pdb, type="launch")
