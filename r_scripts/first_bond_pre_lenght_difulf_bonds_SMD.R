part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
number_y<-c(0:10000)
part<-part_start
setwd(part)
v_start<-c('SMD')
i<-1
j<-1
for (i in 1:length(v_start)) {
  df_structure<-read.csv("df_structure.csv",stringsAsFactors = F)
#  df_structure<-df_structure%>%filter(type=="extra")
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
  if(!dir.exists(paste0(v_start[i],"/din/lenght"))){dir.create(paste0(v_start[i],"/din/lenght"))}
  df_structure<-df_structure%>%select(type,colnames(df_structure)[c(ncol(df_structure)-2,ncol(df_structure)-1)],bonds)
  colnames(df_structure)<-c("type", "C1", "C2", "bonds")
  for (j in 1:nrow(df_structure)) { 
    df_bonds<-df_structure%>%filter(C1==df_structure$C1[j])
    df_bonds<-df_bonds%>%filter(C2==df_structure$C2[j])
    for (q in 1:2) {
      if(!file.exists(paste0(v_start[i],"/din/lenght/lenght_bonds_",df_structure$bonds[j],"_",q,".csv"))){
        frame_number<-length(list.files(path = paste0(v_start[i],"/din/pdb_second/",df_structure$bonds[j],"_",q)))-1
        df_bonds<-df_bonds%>%mutate(c=NA)
        df_bonds<-df_bonds%>%mutate(bonds=df_structure$bonds[j])
        df_bonds<-unique(df_bonds)
        if(frame_number>0){
          df_lenght<-data.frame(matrix(ncol = 3,nrow = (frame_number+1)))
          colnames(df_lenght)<-c("frame","lenght","c")
          df_lenght$frame<-c(0:frame_number)
          df_lenght<-left_join(df_lenght,df_bonds,by="c")
          df_lenght$c<-NULL
          for (p in 1:nrow(df_lenght)) {
            pdb<-read.pdb(paste0(v_start[i],"/din/pdb_second/",df_structure$bonds[j],"_",q,"/frame_",df_lenght$frame[p],".pdb"))
            df_pdbs<-pdb$atom
            df_pdb_C1<-df_pdbs%>%filter(resno==df_lenght$C1[p])
            df_pdb_C2<-df_pdbs%>%filter(resno==df_lenght$C2[p])
            df_lenght$lenght[p]<-sqrt((df_pdb_C1$x[1]-df_pdb_C2$x[1])^2+(df_pdb_C1$y[1]-df_pdb_C2$y[1])^2+(df_pdb_C1$z[1]-df_pdb_C2$z[1])^2)
          }
          df_lenght<-df_lenght%>%mutate(attempt=q)
          df_lenght<-df_lenght%>%mutate(method=v_start[i])
          df_lenght<-df_lenght%>%mutate(system=df_structure$bonds[j])
          write.csv(df_lenght,file = paste0(v_start[i],"/din/lenght/lenght_bonds_",df_structure$bonds[j],"_",q,".csv"),row.names = F)
        }
      }
    }
  }
}
