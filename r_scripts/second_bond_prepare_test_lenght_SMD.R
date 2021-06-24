#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)
library(bio3d)
setwd(part_start)
part<-paste0(part_start,'second_bond/')
setwd(part)
df_structure<-read.csv(paste0(part,'df_structure.csv'),stringsAsFactors = F)
number_y<-c(0:10000)

v_start<-c("SMD","SMD_rev")
i<-1
j<-1
q<-1
p<-1
find_length_bonds<-function(df_structure,v_method,v_system,v_attempt){
  df_leng_test<-df_structure%>%mutate(min_lenght=NA)
  df_leng_test<-df_leng_test%>%mutate(attempt=v_attempt)
  df_leng_test<-df_leng_test%>%mutate(method=v_method)
  df_lenght<-read.csv( paste0(v_method,"/din/lenght/lenght_bonds_loop_",v_system,"_",v_attempt,".csv"),stringsAsFactors = F)
  df_lenght<-df_lenght%>%mutate(attempt=v_attempt)
  df_lenght<-df_lenght%>%mutate(method=v_method)
  df_lenght<-df_lenght%>%mutate(system=v_system)
  df_lenght<-df_lenght%>%filter(first_bond<4)
  df_lenght<-df_lenght%>%filter(second_bond<6)
  df_leng_test<-right_join(df_leng_test,df_lenght,by = c("attempt", "method", "bonds"="system"))
  return(df_leng_test)
}
for (i in 1:length(v_start)) {
  if(!dir.exists(v_start[i])){dir.create(v_start[i])}
  if(!dir.exists(paste0(v_start[i],"/protein"))){dir.create(paste0(v_start[i],"/protein"))}
  if(!dir.exists(paste0(v_start[i],"/protein/tcl"))){dir.create(paste0(v_start[i],"/protein/tcl"))}
  if(!dir.exists(paste0(v_start[i],"/pdb"))){dir.create(paste0(v_start[i],"/pdb"))}
  if(!dir.exists(paste0(v_start[i],"/quench"))){dir.create(paste0(v_start[i],"/quench"))}
  if(!dir.exists(paste0(v_start[i],"/din"))){dir.create(paste0(v_start[i],"/din"))}
  if(!dir.exists(paste0(v_start[i],"/din/lenght"))){dir.create(paste0(v_start[i],"/din/lenght"))}
  if(!dir.exists(paste0(v_start[i],"/din/min_lenght"))){dir.create(paste0(v_start[i],"/din/min_lenght"))}
  for (j in 1:nrow(df_structure)) {
    v_method<-v_start[i]
    v_system<-df_structure$bonds[j]
    df_length_bonds_1<-find_length_bonds(df_structure,v_method,v_system,v_attempt=1)
    df_length_bonds_2<-find_length_bonds(df_structure,v_method,v_system,v_attempt=2)
    df_length_bonds<-rbind(df_length_bonds_1,df_length_bonds_2)
    write.csv(df_length_bonds,paste0(v_start[i],"/din/min_lenght/lenght_",v_system,".csv"),row.names = F)
  }
  df_leng_test<-read.csv(paste0(v_start[i],"/din/min_lenght/lenght_",df_structure$bonds[1],".csv"),stringsAsFactors =  F)
  for (j in 2:nrow(df_structure)) {
    df_leng_test_add<-read.csv(paste0(v_start[i],"/din/min_lenght/lenght_",df_structure$bonds[j],".csv"),stringsAsFactors = F)
    df_leng_test<-rbind(df_leng_test,df_leng_test_add)
  }
  write.csv(df_leng_test,paste0(part_start,'second_bond/',v_start[i],"_length.csv"),row.names = F)
}
df_leng_test<-read.csv(paste0(part_start,'second_bond/',v_start[1],"_length.csv"))
df_leng_test_add<-read.csv(paste0(part_start,'second_bond/',v_start[2],"_length.csv"))
df_leng_test<-rbind(df_leng_test,df_leng_test_add)
df_leng_test<-df_leng_test
df_leng_test<-df_leng_test%>%filter(min_lenght<6)
df_leng_test<-df_leng_test%>%group_by(bonds)%>%mutate(lenght=min(min_lenght))%>%filter(lenght==min_lenght)
write.csv(df_leng_test,paste0(part_start,"second_bond/min_length.csv"),row.names = F)
if(!dir.exists(paste0(part,"protein_SMD"))){dir.create(paste0(part,"protein_SMD"))}
i<-1
for (i in 1:nrow(df_leng_test)){
  pdb<-read.pdb(paste0(part,df_leng_test$method[i],"/din/pdb_second/NaPi2b_",df_leng_test$system[i],"_",df_leng_test$attempt[i],"/frame_",df_leng_test$frame[i],".pdb"))
  write.pdb(pdb,paste0(part,"protein_SMD/",df_leng_test$bonds[i],".pdb"))
}

