part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
number_y<-c(0:10000)
if (!dir.exists("second_bond")){dir.create("second_bond")}
df_structure<-read.csv("first_bond/df_structure.csv",stringsAsFactors = F)

v_start<-c('SMD')
#df_structure<-df_structure%>%mutate(first_bond=paste0(C1,"-",C2))
for (i in 1:nrow(df_structure)) {
  if(!file.exists(paste0(part_start,"first_bond/stablilisation/pdb/",df_structure$bonds[i],".coor"))){
    df_structure$bonds[i]<-NA
  }
}
df_structure<-df_structure%>%filter(!is.na(bonds))
df_structure_add<-read.csv("first_bond/df_structure.csv",stringsAsFactors = F)

df_structure<-df_structure%>%mutate(type="glue")
df_structure_add<-df_structure_add%>%mutate(type="glue")

df_structure$bonds<-NULL
df_structure_add$bonds<-NULL

df_structure<-full_join(df_structure,df_structure_add,by="type")
colnames(df_structure)<-c("type",paste0("C",c(1:(ncol(df_structure)-1))))
df_structure<-df_structure%>%filter(abs(C3-C1)>5)
df_structure<-df_structure%>%filter(abs(C3-C2)>5)
df_structure<-df_structure%>%filter(abs(C4-C1)>5)
df_structure<-df_structure%>%filter(abs(C4-C2)>5)
#df_structure<-rbind(df_structure_1,df_structure_2)
df_structure<-df_structure%>%mutate(bonds=paste0(C1,"-",C2,"_",C3,"-",C4))
write.csv(df_structure,paste0("second_bond/df_structure.csv"),row.names = F)
df_structure<-df_structure%>%mutate(first_bond=paste0(C1,"-",C2))
part<-paste0(part_start,"second_bond/")
if (!dir.exists(part)){dir.create(part)}
if (!dir.exists(paste0(part_start,"second_bond/"))){dir.create(paste0(part_start,"second_bond/"))}
if (!dir.exists(paste0(part_start,"second_bond/protein"))){dir.create(paste0(part_start,"second_bond/protein"))}
if (!dir.exists(paste0(part_start,"first_bond/fin_protein/"))){dir.create(paste0(part_start,"first_bond/fin_protein/"))}
v_structure<-unique(df_structure$first_bond)
for (i in 1:length(v_structure)) {
  pdb<-read.pdb(paste0(part_start,"first_bond/stablilisation/pdb/",v_structure[i],".coor"))
  write.pdb(pdb,paste0(part_start,"first_bond/fin_protein/",v_structure[i],".pdb"))
  write.pdb(pdb,paste0(part_start,"second_bond/protein/",v_structure[i],".pdb"))
}
