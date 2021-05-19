part_start <- commandArgs(trailingOnly=TRUE)
setwd(part_start)
seqv<-read.fasta(paste0("SLC34.fasta"))
seq<-seqaln(seqv, id=NULL, profile=NULL, exefile="muscle", outfile="aln.fa")
v_seq<-seq$id
df_seq<-data.frame(matrix(ncol=4,nrow = length(v_seq)))
colnames(df_seq)<-c("alig_id","seq_id","protein_name","spesies")
df_seq$alig_id<-v_seq
i<-1
for (i in 1:nrow(df_seq)) {
  a<-strsplit(df_seq$alig_id[i],split = "|",fixed = T)[[1]]
  df_seq$seq_id[i]<-a[2]
  b<-strsplit(a[3],split = "_",fixed = T)[[1]]
  df_seq$protein_name[i]<-b[1]
  df_seq$spesies[i]<-b[2]
}
rm(a,b)
df_seq<-df_seq%>%mutate(seq=paste0(protein_name,"_",spesies))
v_ali<-seq$ali
v_ali<-t(v_ali)
df_aligment<-data.frame(v_ali)
colnames(df_aligment)<-df_seq$seq
df_aligment <- data.frame(lapply(df_aligment, as.character), stringsAsFactors=FALSE)
df_aligment<-df_aligment%>%filter(NPT2B_HUMAN!="-")
df_aligment<-df_aligment%>%mutate(resno=1:nrow(df_aligment))
df_aligment<-df_aligment%>%mutate(persent_SLC34=0)
df_aligment<-df_aligment%>%mutate(persent_NaPi2b=0)
df_seq<-df_seq%>%mutate(SLC34_num=1:nrow(df_seq))
df_seq<-df_seq%>%mutate(NaPi2b_num=1:nrow(df_seq))
df_seq$NaPi2b_num[df_seq$protein_name!="NPT2B"]<-NA
df_aligment_NaPi2b<-df_aligment%>%select(NPT2B_HUMAN, NPT2B_RAT, NPT2B_MOUSE, NPT2B_BOVIN, NPT2B_PONAB,number,persent_SLC34,persent_NaPi2b)
for (i in 1:(nrow(df_seq))) {
  df_aligment$persent_SLC34[df_aligment$NPT2B_HUMAN==df_aligment[,i]]<-df_aligment$persent_SLC34[df_aligment$NPT2B_HUMAN==df_aligment[,i]]+1
  if (!is.na(df_seq$NaPi2b_num[i])) {
    df_aligment$persent_NaPi2b[df_aligment$NPT2B_HUMAN==df_aligment[,i]]<-df_aligment$persent_NaPi2b[df_aligment$NPT2B_HUMAN==df_aligment[,i]]+1
  }
}
#df_seq[!is.na(df_seq$NaPi2b_num),]
v_NaPi2b<-(nrow(df_seq[!is.na(df_seq$NaPi2b_num),]))
v_SLC34<-(nrow(df_seq))
df_aligment<-df_aligment%>%mutate(persent_SLC34=persent_SLC34/v_SLC34*100)
df_aligment<-df_aligment%>%mutate(persent_NaPi2b=persent_NaPi2b/v_NaPi2b*100)
df_aligment<-df_aligment%>%select(resno, persent_SLC34,persent_NaPi2b)
write.csv(df_aligment,"fin_aligment.csv",row.names = F)