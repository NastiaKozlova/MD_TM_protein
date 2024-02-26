part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
v_name<-list.files("MD")
name<-v_name[1]
v_namd<-1000
i<-1
for (name in v_name) {
    df_tcl<-read.table(file = paste0("MD/",name,"/namd/step7_production.inp"),sep="+")
    for (i in 1:nrow(df_tcl)) {
        if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = "outputName", fixed = TRUE))){
            finish_file<-i
        }
        if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = "set inputname", fixed = TRUE))){
            set_inputname<-i
        }
        if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = "useFlexibleCell", fixed = TRUE))){
            useFlexibleCell<-i
        }
        if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = "wrapAll", fixed = TRUE))){
            wrapAll<-i
        }
        
        if(length(df_tcl[i,1])>0&(grepl( df_tcl[i,1], pattern = 'hexa', fixed = TRUE))){
            v_hexa<-i
        }
        
    }

    df_conf<-df_tcl
    i<-0
    df_conf[set_inputname,1]<-paste0("set inputname           step7_production;")
    df_conf[useFlexibleCell,1]<-paste0("useFlexibleCell         yes;")
    df_conf[wrapAll,1]<-paste0("wrapAll         no;")
    df_conf[v_hexa,1]<-paste0('if { $boxtype == "hexa" } {')
    df_conf[finish_file,1]<-paste0("outputName              step7.",i+1,"_production; # base name for output from this run\n")
    write.table(df_conf,paste0("MD/",name,"/namd/step7.",i+1,"_production.inp"),row.names=F,col.names=F,quote = F)
    for (j in 1:v_namd) {
#        df_conf<-df_tcl
        df_conf[finish_file,1]<-paste0("outputName              step7.",j+1,"_production; # base name for output from this run\n")
        df_conf[set_inputname,1]<-paste0("set inputname           step7.",j,"_production;")
        write.table(df_conf,paste0("MD/",name,"/namd/step7.",j+1,"_production.inp"),row.names=F,col.names=F,quote = F)
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
