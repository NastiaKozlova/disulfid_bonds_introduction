part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
number_y<-c(0:10000)
#part<-paste0(part_start,'first_bond/')
part<-part_start
setwd(part)
v_start<-c('SMD')
#v_start<-c('SMD')
i<-1
for (i in 1:length(v_start)) {
  df_structure<-read.csv("df_structure.csv",stringsAsFactors = F)
  q<-1
  for (j in 1:nrow(df_structure)) {
    if (!file.exists(paste0(v_start[i],"/quench/quench_",df_structure$bonds[j],"_",q,".dcd"))) {
      df_structure$bonds[j]<-NA
    }
  }
  df_structure<-df_structure%>%filter(!is.na(bonds))
  if(!dir.exists(v_start[i])){dir.create(v_start[i])}
  if(!dir.exists(paste0(v_start[i],"/protein"))){dir.create(paste0(v_start[i],"/protein"))}
  if(!dir.exists(paste0(v_start[i],"/protein/tcl"))){dir.create(paste0(v_start[i],"/protein/tcl"))}
  if(!dir.exists(paste0(v_start[i],"/pdb"))){dir.create(paste0(v_start[i],"/pdb"))}
  if(!dir.exists(paste0(v_start[i],"/quench"))){dir.create(paste0(v_start[i],"/quench"))}
  if(!dir.exists(paste0(v_start[i],"/din"))){dir.create(paste0(v_start[i],"/din"))}
  if(!dir.exists(paste0(v_start[i],"/din/tcl"))){dir.create(paste0(v_start[i],"/din/tcl"))}
  if(!dir.exists(paste0(v_start[i],"/din/pdb_second"))){dir.create(paste0(v_start[i],"/din/pdb_second"))}
  for (j in 1:nrow(df_structure)) {
    for (q in 1:2) {
      if(file.exists(paste0(v_start[i],"/quench/quench_",df_structure$bonds[j],"_",q,".dcd"))){
      if (!dir.exists(paste0(v_start[i],"/din/pdb_second/",df_structure$bonds[j],"_",q))) {dir.create(paste0(v_start[i],"/din/pdb_second/",df_structure$bonds[j],"_",q))}
      df_tcl<-data.frame(matrix(nrow = 1,ncol = 7))
      df_tcl[1,1]<-paste0("cd ", part,v_start[i],"\nmol new {protein/",df_structure$bonds[j],".psf} type {psf}")
      df_tcl[1,2]<-paste0("mol addfile {quench/quench_",df_structure$bonds[j],"_",q,".dcd} type {dcd} first 0 last -1 step 1 waitfor all")
      df_tcl[1,3]<-paste0("set nf [molinfo top get numframes]")
      df_tcl[1,4]<-paste0("for {set i 0 } {$i < $nf} {incr i} {")
      df_tcl[1,5]<-paste0('[atomselect top "protein and name SG" frame $i] writepdb din/pdb_second/',df_structure$bonds[j],'_',q,'/frame_$i.pdb')
      df_tcl[1,6]<-paste0("}")
      df_tcl[1,7]<-"mol delete all"
      df_tcl[1,8]<-"exit now"
      write.table(df_tcl,paste0(v_start[i],"/din/tcl/Second_str_",df_structure$bonds[j],"_",q,".tcl"),sep="\n",col.names = F,row.names = F,quote = F)
      system(command = paste0("vmd -dispdev text -e ", v_start[i],"/din/tcl/Second_str_",df_structure$bonds[j],"_",q,".tcl"),ignore.stdout=T,wait = T)
    }
  }
}
}
