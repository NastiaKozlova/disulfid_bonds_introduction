part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
number_y<-c(0:10000)
part<-part_start
setwd(part)
v_start<-c('SMD')

df_structure<-read.csv("df_structure.csv",stringsAsFactors = F)
if(!dir.exists(paste0(part,"tcl_pdb_find/"))){dir.create(paste0(part,"tcl_pdb_find/"))}
if(!dir.exists(paste0(part,"protein_SMD/"))){dir.create(paste0(part,"protein_SMD/"))}

i<-1
j<-1
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
  if(!dir.exists(paste0(v_start[i],"/din/min_lenght"))){dir.create(paste0(v_start[i],"/din/min_lenght"))}
  if(nrow(df_structure)>1){
    for (j in 1:nrow(df_structure)) {
      df_lenght_1<-read.csv(paste0(v_start[i],"/din/lenght/lenght_bonds_",df_structure$bonds[j],"_",1,".csv"),stringsAsFactors = F)
      df_lenght_2<-read.csv(paste0(v_start[i],"/din/lenght/lenght_bonds_",df_structure$bonds[j],"_",2,".csv"),stringsAsFactors = F)
      df_lenght_1<-df_lenght_1%>%filter(lenght<5)
      df_lenght_2<-df_lenght_2%>%filter(lenght<5)
      df_lenght<-df_lenght_2
      if(nrow(df_lenght_1)>nrow(df_lenght_2)){df_lenght<-df_lenght_1}
      if(nrow(df_lenght)>0){
        write.csv(df_lenght,file = paste0(v_start[i],"/din/min_lenght/",df_structure$bonds[j],".csv"),row.names = F)
      }
    }
  }
}
i<-2
for (i in 1:length(v_start)) {
  df_structure<-read.csv("df_structure.csv",stringsAsFactors = F)
  if(!dir.exists(paste0(v_start[i],"/din/tcl_pdb_find"))){dir.create(paste0(v_start[i],"/din/tcl_pdb_find"))}
  df_structure<-df_structure%>%filter(!is.na(bonds))
  for (j in 1:nrow(df_structure)) {
    
    if (!file.exists(paste0(v_start[i],"/din/min_lenght/",df_structure$bonds[j],".csv"))) {
      df_structure$bonds[j]<-NA
    }
  }
  df_structure<-df_structure%>%filter(!is.na(bonds))
  if(nrow(df_structure)>0){
  df_lenght<-read.csv(paste0(v_start[i],"/din/min_lenght/",df_structure$bonds[1],".csv"),stringsAsFactors = F)
  df_lenght<-df_lenght%>%mutate(number=nrow(df_lenght))
  if(nrow(df_structure)>1){
    for (j in 2:nrow(df_structure)) {
      df_lenght_add<-read.csv(paste0(v_start[i],"/din/min_lenght/",df_structure$bonds[j],".csv"),stringsAsFactors = F)
      df_lenght_add<-df_lenght_add%>%mutate(number=nrow(df_lenght_add))
      df_lenght<-rbind(df_lenght,df_lenght_add)
    }
  }
  df_lenght<-df_lenght%>%group_by(bonds)%>%mutate(min_frame=min(frame))
  df_lenght<-ungroup(df_lenght)
  df_lenght<-df_lenght%>%filter(min_frame==frame)
  write.csv(df_lenght,paste0(v_start[i],"_min_length.csv"),row.names = F)
  }
}
df_lenght<-read.csv(paste0(v_start[1],"_min_length.csv"))
if(length(v_start)>1){
  for (i in 2:length(v_start)) {
    if(file.exists(paste0(v_start[i],"_min_length.csv"))){
      df_lenght_add<-read.csv(paste0(v_start[i],"_min_length.csv"),stringsAsFactors = F)
      df_lenght<-rbind(df_lenght,df_lenght_add)
    }
  }
}
df_lenght<-df_lenght%>%group_by(bonds)%>%mutate(min_number=min(number))
df_lenght<-ungroup(df_lenght)
df_lenght<-df_lenght%>%filter(min_number==number)
write.csv(df_lenght,"min_length.csv",row.names = F)

for (j in 1:nrow(df_lenght)) {
  df_tcl<-data.frame(matrix(nrow = 1,ncol = 5))
  df_tcl[1,1]<-paste0("cd ", part,df_lenght$method[j],"\nmol new {protein/",df_lenght$bonds[j],".psf} type {psf}")
  df_tcl[1,2]<-paste0("mol addfile {quench/quench_",df_lenght$bonds[j],"_",df_lenght$attempt[1],".dcd} type {dcd} first 0 last -1 step 1 waitfor all")
  df_tcl[1,3]<-paste0('[atomselect top all frame ',df_lenght$frame[j],'] ', 'writepdb ',part,'protein_SMD/',df_lenght$bonds[j],'.pdb')
  df_tcl[1,4]<-"mol delete all\n\n"
  df_tcl[1,5]<-"exit now"
  write.table(df_tcl,paste0(part,"tcl_pdb_find/Second_str_",df_lenght$bonds[j],".tcl"),sep="\n",col.names = F,row.names = F,quote = F)
}
for (j in 1:nrow(df_lenght)) {
  print(paste0("vmd -dispdev text -e ",part,"tcl_pdb_find/Second_str_",df_lenght$bonds[j],".tcl"))
  system(command =paste0("vmd -dispdev text -e ",part,"tcl_pdb_find/Second_str_",df_lenght$bonds[j],".tcl"),ignore.stdout=T,wait = T)
}
