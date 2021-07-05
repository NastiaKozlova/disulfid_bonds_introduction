part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
number_y<-c(0:10000)
if (!dir.exists("third_bond")){dir.create("third_bond")}
if (!dir.exists("analysis")){dir.create("analysis")}

#find most common first bond
df_structure<-read.csv("second_bond/df_structure.csv",stringsAsFactors = F)
df_structure<-df_structure%>%mutate(number=1)
for (i in 1:nrow(df_structure)) {
  if(!file.exists(paste0(part_start,"second_bond/stablilisation/pdb/",df_structure$bonds[i],".coor"))){
    df_structure$bonds[i]<-NA
    df_structure$number[i]<-0
  }
}
df_structure_second<-df_structure



df_structure<-read.csv("first_bond/df_structure.csv",stringsAsFactors = F)
df_structure<-df_structure%>%mutate(number=1)
for (i in 1:nrow(df_structure)) {
  if(!file.exists(paste0(part_start,"first_bond/stablilisation/pdb/",df_structure$bonds[i],".coor"))){
    df_structure$bonds[i]<-NA
    df_structure$number[i]<-0
  }
}
df_structure_first<-df_structure
df_structure_first<-df_structure_first%>%mutate(C3=C1)
df_structure_first<-df_structure_first%>%mutate(C4=C2)
df_structure_first<-df_structure_first%>%mutate(C1=0)
df_structure_first<-df_structure_first%>%mutate(C2=0)
df_structure_first<-df_structure_first%>%select(type,C1,C2,C3,C4,bonds,number)
df_structure<-rbind(df_structure_first,df_structure_second)
df_structure<-df_structure%>%mutate(first_bond=paste0(C1,"-",C2))
df_structure<-df_structure%>%mutate(second_bond=paste0(C3,"-",C4))
df_structure_bonds<-df_structure%>%filter(number==1)
p<-ggplot(data=df_structure_bonds)+labs(x = "First bond",y="Second bond")+
  geom_point(aes(y=first_bond,x=second_bond,colour=number))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=-90),panel.grid.major.x = element_line(colour = "black"))+ 
  guides(colour = "none")

v_heigth<-length(unique(df_structure_bonds$first_bond))
v_whith<-length(unique(df_structure_bonds$second_bond))
ggsave(p,filename = paste0("analysis/df_first_second_bond.png"), width = (v_whith)*0.5, height =v_heigth*0.5, units = c("cm"), dpi = 200 )
df_structure<-df_structure%>%group_by(first_bond)%>%mutate(start_bonds=n())

df_structure<-df_structure%>%mutate(finish_bonds=0)
df_structure_test<-df_structure%>%filter(number==1)
df_structure_test<-df_structure_test%>%group_by(first_bond)%>%mutate(finish_bonds=n())
i<-1
for (i in 1:nrow(df_structure_test)) {
  df_structure$finish_bonds[df_structure$bonds%in%df_structure_test$bonds[i]]<-df_structure_test$finish_bonds[i]
}
df_structure<-df_structure%>%mutate(persent_bonds=finish_bonds/start_bonds*100)
df_structure<-ungroup(df_structure)
df_structure<-df_structure%>%group_by(first_bond)%>%mutate(max_persent=max(persent_bonds))
df_structure_a<-df_structure%>%filter(persent_bonds==max_persent)
df_structure_a<-df_structure_a%>%select(first_bond,persent_bonds)
df_structure_first_bond<-unique(df_structure_a)
df_structure_first_bond<-df_structure_first_bond%>%filter(persent_bonds>0)
write.csv(df_structure_first_bond,paste0("analysis/df_structure_first_bond.csv"),row.names = F)
p<-ggplot(data=df_structure_first_bond)+
  labs(y="First bonds",x="Persent interactions")+
  geom_point(aes(y=first_bond,x=persent_bonds))+
  theme_bw()+
  theme(axis.text.x=element_text(angle=-90),panel.grid.major.x = element_line(colour = "black"))
v_heigth<-length(unique(df_structure_first_bond$first_bond))
#v_whith<-length(unique(df_structure_first_bond$second_bond))
ggsave(p,filename = paste0("analysis/df_structure_first_bond.png"), width = 20, height =v_heigth*0.5, units = c("cm"), dpi = 200 )
#find most common second bond
df_structure<-read.csv("second_bond/df_structure.csv",stringsAsFactors = F)
df_structure<-df_structure%>%mutate(number=1)
for (i in 1:nrow(df_structure)) {
  if(!file.exists(paste0(part_start,"second_bond/stablilisation/pdb/",df_structure$bonds[i],".coor"))){
    df_structure$bonds[i]<-NA
    df_structure$number[i]<-0
  }
}
df_structure_second<-df_structure



df_structure<-read.csv("first_bond/df_structure.csv",stringsAsFactors = F)
df_structure<-df_structure%>%mutate(number=1)
for (i in 1:nrow(df_structure)) {
  if(!file.exists(paste0(part_start,"first_bond/stablilisation/pdb/",df_structure$bonds[i],".coor"))){
    df_structure$bonds[i]<-NA
    df_structure$number[i]<-0
  }
}
df_structure_first<-df_structure
df_structure_first<-df_structure_first%>%mutate(C3=C1)
df_structure_first<-df_structure_first%>%mutate(C4=C2)
df_structure_first<-df_structure_first%>%mutate(C1=0)
df_structure_first<-df_structure_first%>%mutate(C2=0)
df_structure_first<-df_structure_first%>%select(type,C1,C2,C3,C4,bonds,number)
df_structure<-rbind(df_structure_first,df_structure_second)
#df_structure<-df_structure%>%filter(number==1)
df_structure<-df_structure%>%mutate(first_bond=paste0(C1,"-",C2))
df_structure<-df_structure%>%mutate(second_bond=paste0(C3,"-",C4))
df_structure_bonds<-df_structure%>%filter(number==1)
ggplot(data=df_structure_bonds)+geom_point(aes(y=first_bond,x=second_bond,colour=number))+
  theme(axis.text.x=element_text(angle=-90))

df_structure<-df_structure%>%group_by(second_bond)%>%mutate(start_bonds=n())

df_structure<-df_structure%>%mutate(finish_bonds=0)
df_structure_test<-df_structure%>%filter(number==1)
df_structure_test<-df_structure_test%>%group_by(second_bond)%>%mutate(finish_bonds=n())
i<-1
for (i in 1:nrow(df_structure_test)) {
  df_structure$finish_bonds[df_structure$bonds%in%df_structure_test$bonds[i]]<-df_structure_test$finish_bonds[i]
}
df_structure<-df_structure%>%mutate(persent_bonds=finish_bonds/start_bonds*100)
df_structure<-ungroup(df_structure)
df_structure<-df_structure%>%group_by(second_bond)%>%mutate(max_persent=max(persent_bonds))
df_structure_a<-df_structure%>%filter(persent_bonds==max_persent)
df_structure_a<-df_structure_a%>%select(second_bond,persent_bonds)
df_structure_a<-unique(df_structure_a)
df_structure_second_bond<-df_structure_a%>%filter(persent_bonds>0)
write.csv(df_structure_first_bond,paste0("analysis/df_structure_second_bond.csv"),row.names = F)
#p<-ggplot(data=df_structure_second_bond)+geom_point(aes(y=second_bond,x=persent_bonds))+
#  theme(axis.text.x=element_text(angle=-90))
#ggsave(p,filename = paste0("analysis/df_structure_second_bond.png"))
p<-ggplot(data=df_structure_second_bond)+
  geom_point(aes(y=second_bond,x=persent_bonds))+
  labs(y="Second bonds",x="Persent interactions")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=-90),panel.grid.major.x = element_line(colour = "black"))
#v_heigth<-length(unique(df_structure_first_bond$first_bond))
v_whith<-length(unique(df_structure_second_bond$second_bond))
ggsave(p,filename = paste0("analysis/df_structure_second_bond.png"), height= 20, width =v_heigth*0.5, units = c("cm"), dpi = 200 )






#find most common combinations between first and second bonds
df_structure<-read.csv("second_bond/df_structure.csv",stringsAsFactors = F)
for (i in 1:nrow(df_structure)) {
  if(!file.exists(paste0(part_start,"second_bond/stablilisation/pdb/",df_structure$bonds[i],".coor"))){
    df_structure$bonds[i]<-NA
  }
}
df_structure<-df_structure%>%filter(!is.na(bonds))
df_structure<-df_structure%>%mutate(first_bond=paste0(C1,"-",C2))
df_structure<-df_structure%>%mutate(second_bond=paste0(C3,"-",C4))

df_structure_first_bond<-df_structure_first_bond%>%filter(persent_bonds>quantile(df_structure_first_bond$persent_bonds,probs = 0.75))
df_structure_second_bond<-df_structure_second_bond%>%filter(persent_bonds>quantile(df_structure_second_bond$persent_bonds,probs = 0.5))

df_structure<-df_structure[df_structure$first_bond%in%df_structure_first_bond$first_bond,]
df_structure<-df_structure[df_structure$second_bond%in%df_structure_second_bond$second_bond,]
write.csv(df_structure,"analysis/df_structure_fin.csv",row.names = F)
#df_structure$bonds<-NULL
df_structure_add<-df_structure%>%select(type,C3,C4,C1,C2,bonds)
df_structure<-df_structure%>%select(type,C1,C2,C3,C4,bonds)
colnames(df_structure_add)<-colnames(df_structure)
df_structure<-rbind(df_structure,df_structure_add)
#df_structure<-df_structure%>%mutate(second_bond=paste0(C1,"-",C2))

df_structure_mod<-df_structure
#df_structure_mod<-df_structure%>%filter(C1<C3)
for (i in 1:nrow(df_structure)) {
  if(df_structure$C1[i]>df_structure$C3[i]){
    df_structure_mod$C1[i]<-df_structure$C3[i]
    df_structure_mod$C2[i]<-df_structure$C4[i]
    df_structure_mod$C3[i]<-df_structure$C1[i]
    df_structure_mod$C4[i]<-df_structure$C2[i]
  }
}
df_structure<-df_structure_mod%>%mutate(bonds=paste0(C1,"-",C2,"_",C3,"-",C4))
df_structure<-unique(df_structure)
write.csv(df_structure,"second_bond/df_structure_fin.csv",row.names = F)

df_structure<-df_structure%>%mutate(first_bond=paste0(C1,"-",C2))
df_structure<-df_structure%>%mutate(second_bond=paste0(C3,"-",C4))
df_structure<-df_structure%>%arrange(first_bond)
df_structure<-unique(df_structure)
df_structure<-ungroup(df_structure)
df_structure<-df_structure%>%select(type,C1,C2,C3,C4)
df_structure_add<-read.csv("first_bond/df_structure.csv",stringsAsFactors = F)

df_structure<-df_structure%>%mutate(type="glue")
df_structure_add<-df_structure_add%>%mutate(type="glue")

df_structure$bonds<-NULL
df_structure_add$bonds<-NULL
df_structure$number<-NULL

df_structure<-full_join(df_structure,df_structure_add,by="type")
colnames(df_structure)<-c("type",paste0("C",c(1:(ncol(df_structure)-1))))
df_structure<-df_structure%>%filter(abs(C3-C1)>5)
df_structure<-df_structure%>%filter(abs(C3-C2)>5)
df_structure<-df_structure%>%filter(abs(C4-C1)>5)
df_structure<-df_structure%>%filter(abs(C4-C2)>5)
df_structure<-df_structure%>%filter(abs(C5-C1)>5)
df_structure<-df_structure%>%filter(abs(C5-C2)>5)
df_structure<-df_structure%>%filter(abs(C5-C3)>5)
df_structure<-df_structure%>%filter(abs(C5-C4)>5)
df_structure<-df_structure%>%filter(abs(C6-C1)>5)
df_structure<-df_structure%>%filter(abs(C6-C2)>5)
df_structure<-df_structure%>%filter(abs(C6-C3)>5)
df_structure<-df_structure%>%filter(abs(C6-C4)>5)
df_structure<-df_structure%>%mutate(bonds=paste0(C1,"-",C2,"_",C3,"-",C4,"_",C5,"-",C6))
write.csv(df_structure,paste0("third_bond/second_bonds_filtered.csv"),row.names = F)

df_structure_add<-left_join(df_structure,df_structure_mod,by=c("type", "C1", "C2", "C3", "C4"))
df_structure_add<-df_structure_add%>%select(type,C1,C2,C3,C4,bonds.y)
df_structure_add<-unique(df_structure_add)
df_structure_add<-df_structure_add%>%mutate(first_bond=paste0(C1,"-",C2,"_",C3,"-",C4))
df_structure_add<-df_structure_add%>%group_by(first_bond)%>%mutate(number=n())
df_structure_1<-df_structure_add%>%filter(number==1)
df_structure_2<-df_structure_add%>%filter(number==2)
df_structure_2<-df_structure_2%>%filter(first_bond==bonds.y)
df_structure_test<-rbind(df_structure_1,df_structure_2)
#df_structure_TEMP<-left_join(df_structure,df_structure_test,by=c("type", "C1", "C2", "C3", "C4"))

part<-paste0(part_start,"third_bond/")
if (!dir.exists(part)){dir.create(part)}
if (!dir.exists(paste0(part_start,"third_bond/"))){dir.create(paste0(part_start,"third_bond/"))}
if (!dir.exists(paste0(part_start,"third_bond/protein"))){dir.create(paste0(part_start,"third_bond/protein"))}
if (!dir.exists(paste0(part_start,"second_bond/fin_protein/"))){dir.create(paste0(part_start,"second_bond/fin_protein/"))}
for (i in 1:nrow(df_structure_test)) {
  pdb<-read.pdb(paste0(part_start,"second_bond/stablilisation/pdb/",df_structure_test$bonds.y[i],".coor"))
  write.pdb(pdb,paste0(part_start,"second_bond/fin_protein/",df_structure_test$first_bond[i],".pdb"))
  write.pdb(pdb,paste0(part_start,"third_bond/protein/",df_structure_test$first_bond[i],".pdb"))
}

