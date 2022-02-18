part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
v_name<-list.files("MD")
name<-v_name[1]
v_namd<-1000
for (name in v_name) {
  df_tcl<-read.table(file = paste0("MD/",name,"/namd/step7_production.inp"),sep="?")
  for (i in 1:nrow(df_tcl)) {
    if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = "outputName", fixed = TRUE))){
      finish_file<-i
    }
    if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = "set inputname", fixed = TRUE))){
      set_inputname<-i
    }
  }
  #  grepl( df_tcl[1,i], "outputName", fixed = TRUE)
  
  df_conf<-df_tcl
  i<-0
  df_conf[set_inputname,1]<-paste0("set inputname           step7_production;\n")
  df_conf[finish_file,1]<-paste0("outputName              step7.",i+1,"_production; # base name for output from this run\n")
  write.table(df_conf,paste0("MD/",name,"/namd/step7.",i+1,"_production.inp"),row.names=F,col.names=F,quote = F)
  for (i in 1:v_namd) {
    df_conf<-df_tcl
    df_conf[finish_file,1]<-paste0("outputName              step7.",i+1,"_production; # base name for output from this run\n")
    df_conf[set_inputname,1]<-paste0("set inputname           step7.",i,"_production;\n")
    write.table(df_conf,paste0("MD/",name,"/namd/step7.",i+1,"_production.inp"),row.names=F,col.names=F,quote = F)
  }
}
df_namd<-data.frame(matrix(ncol=2,nrow = 1007))
colnames(df_namd)<-c("comand","conf")
df_namd$comand<-"run_namd"
a<-c( paste0("step6.",c(1:6),"_equilibration.inp > step6.",c(1:6),"_equilibration.out"),
      "step7_production.inp > step7_production.out",
      paste0("step7.",c(1:1000),"_production.inp > step7.",c(1:1000),"_production.out"))
df_namd$conf<-a
df_namd_add<-data.frame(matrix(ncol=2,nrow = 1))
colnames(df_namd_add)<-c("comand","conf")
df_namd_add$comand<-paste0("cd ",part_start,"namd")
df_namd<-rbind(df_namd_add,df_namd)
write.table(df_namd,paste0("run_namd.txt"),sep = " ",row.names = F,col.names = F,quote = F,na = "")
