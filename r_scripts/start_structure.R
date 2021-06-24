part_start = commandArgs(trailingOnly=TRUE)
setwd(part_start)
library(bio3d)
library(dplyr)
num<-5
pdb<-read.pdb(paste0("protein/start.pdb"))
if (!dir.exists("first_bond")){dir.create("first_bond")}
if (!dir.exists("first_bond/protein")){dir.create("first_bond/protein")}
if (!dir.exists("second_bond")){dir.create("second_bond")}
if (!dir.exists("third_bond")){dir.create("third_bond")}
write.pdb(pdb,"first_bond/protein/start.pdb")
df_pdb<-pdb$atom
df_pdb<-df_pdb%>%filter(resid=="CYS")
df_pdb<-df_pdb%>%filter(elety=="CA")
#df_pdb<-df_pdb%>%mutate(name=NA)
df_pdb<-df_pdb%>%mutate(type="TMD")

df_pdb$type[df_pdb$z<(-18)]<-"cyto"
df_pdb$type[df_pdb$z>(18)]<-"extra"
df_pdb_add<-df_pdb%>%filter(resno==350)
df_pdb_add$type<-"extra"
df_pdb<-rbind(df_pdb,df_pdb_add)
df_pdb1<-full_join(df_pdb,df_pdb,by="type")
df_pdb<-df_pdb%>%select(type,resno)
write.csv(df_pdb,"protein/cysteins.csv",row.names=F)
df_pdb<-df_pdb%>%filter(!is.na(resno))
df_pdb1<-full_join(df_pdb,df_pdb,by="type")
colnames(df_pdb1)<-c("type",paste0("C",1:(ncol(df_pdb1)-1)))
df_pdb1<-df_pdb1%>%filter((C1+num)<C2)
df_pdb1<-df_pdb1%>%mutate(bonds=paste0(C1,"-",C2))
write.csv(df_pdb1,paste0("first_bond/df_structure.csv"),row.names = F)
