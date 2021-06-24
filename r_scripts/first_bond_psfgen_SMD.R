part_start = commandArgs(trailingOnly=TRUE)
library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
number_y<-c(0:10000)
part<-paste0(part_start,'first_bond/')
setwd(part)
df_structure<-read.csv("df_structure.csv",stringsAsFactors = F)

v_start<-c('SMD')

i<-2
j<-30
for (i in 1:length(v_start)) {
  if(!dir.exists(v_start[i])){dir.create(v_start[i])}
  if(!dir.exists(paste0(v_start[i],'/protein'))){dir.create(paste0(v_start[i],'/protein'))}
  if(!dir.exists(paste0(v_start[i],'/protein/tcl'))){dir.create(paste0(v_start[i],'/protein/tcl'))}
  if(!dir.exists(paste0(v_start[i],'/pdb'))){dir.create(paste0(v_start[i],'/pdb'))}
  if(!dir.exists(paste0(v_start[i],'/quench'))){dir.create(paste0(v_start[i],'/quench'))}
  if(!dir.exists(paste0(v_start[i],'/din'))){dir.create(paste0(v_start[i],'/din'))}
  for (j in 1:nrow(df_structure)) {
    if (file.exists(paste0('protein/start.pdb'))) {
      fixed_1<-(df_structure$C1[j]-10):(df_structure$C1[j]+10)
      fixed_2<-(df_structure$C2[j]-10):(df_structure$C2[j]+10)
      if(min(fixed_2)-max(fixed_1)<10){
        fixed_1<-(df_structure$C1[j]-10):(df_structure$C2[j]-10)
        fixed_2<-(df_structure$C1[j]+10):(df_structure$C2[j]+10)
      }
      fixed<-unique(c(fixed_1,fixed_2))
      df_tcl<-data.frame(matrix(nrow = 15,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part)
      df_tcl[2,1]<-paste0('mol delete all')
      df_tcl[3,1]<-paste0('package require psfgen')
      df_tcl[4,1]<-paste0('resetpsf')
      df_tcl[5,1]<-paste0('topology ', part_start,'toppar/top_all36_prot.rtf')
      df_tcl[6,1]<-paste0('topology ', part_start,'toppar/toppar_water_ions_namd.str')
      df_tcl[7,1]<-paste0('pdbalias resid ue HIS HSE')
      df_tcl[8,1]<-paste0('pdbalias atom ILE CD1 CD')
      df_tcl[9,1]<-paste0('segment U { pdb protein/start.pdb')
      df_tcl[10,1]<-paste0('}')
      df_tcl[11,1]<-paste0('coordpdb protein/start.pdb U')
      df_tcl[12,1]<-paste0('regenerate angles dihedrals')
      df_tcl[13,1]<-paste0('guesscoord')
      df_tcl[14,1]<-paste0('writepdb ',v_start[i],'/protein/',df_structure$bond[j],'.pdb')
      df_tcl[15,1]<-paste0('writepsf ',v_start[i],'/protein/',df_structure$bond[j],'.psf')
      write.table(df_tcl,paste0(v_start[i],'/protein/tcl/psfgen_',df_structure$bond[j],'.tcl'),row.names = F,col.names = F,quote = F)
      
      df_tcl<-data.frame(matrix(nrow = 10,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part,v_start[i],'/protein')
      df_tcl[2,1]<-paste0('mol delete all')
      df_tcl[3,1]<-paste0('mol new ',df_structure$bonds[j],'.psf')
      df_tcl[4,1]<-paste0('mol addfile ',df_structure$bonds[j],'.pdb')
      df_tcl[5,1]<-paste0('set all [atomselect top "all"]')
      df_tcl[6,1]<-paste0('$all set beta 1')
      df_tcl[7,1]<-paste0('set unfixed [atomselect top " protein and (resid ',paste0(fixed,collapse = " "),')"]')  
      df_tcl[8,1]<-paste0('$unfixed set beta 0')
      df_tcl[9,1]<-paste0('set fixed [atomselect top " protein and (resid  ',df_structure$C1[j],')"]')
      df_tcl[10,1]<-paste0('$fixed set beta 1')
      df_tcl[11,1]<-paste0('$all writepdb ',df_structure$bonds[j],'_1.fix')
      df_tcl[12,1]<-paste0('mol delete all')
      write.table(df_tcl,paste0(v_start[i],'/protein/tcl/fixatom_',df_structure$bond[j],'_1.tcl'),row.names = F,col.names = F,quote = F)
      df_tcl<-data.frame(matrix(nrow = 10,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part,v_start[i],'/protein')
      df_tcl[2,1]<-paste0('mol delete all')
      df_tcl[3,1]<-paste0('mol new ',df_structure$bonds[j],'.psf')
      df_tcl[4,1]<-paste0('mol addfile ',df_structure$bonds[j],'.pdb')
      df_tcl[5,1]<-paste0('set all [atomselect top "all"]')
      df_tcl[6,1]<-paste0('$all set beta 1')
      df_tcl[7,1]<-paste0('set unfixed [atomselect top " protein and (resid ',paste0(fixed,collapse = " "),')"]') 
      df_tcl[8,1]<-paste0('$unfixed set beta 0')
      df_tcl[9,1]<-paste0('set fixed [atomselect top " protein and (resid  ',df_structure$C2[j],')"]')
      df_tcl[10,1]<-paste0('$fixed set beta 1')
      df_tcl[11,1]<-paste0('$all writepdb ',df_structure$bonds[j],'_2.fix')
      df_tcl[12,1]<-paste0('mol delete all')
      write.table(df_tcl,paste0(v_start[i],'/protein/tcl/fixatom_',df_structure$bonds[j],'_2.tcl'),row.names = F,col.names = F,quote = F)
      
      df_tcl<-data.frame(matrix(nrow = 10,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part,v_start[i],'/protein')
      df_tcl[2,1]<-paste0('mol delete all')
      df_tcl[3,1]<-paste0('mol new ',df_structure$bond[j],'.psf')
      df_tcl[4,1]<-paste0('mol addfile ',df_structure$bond[j],'.pdb')
      df_tcl[5,1]<-paste0('set all [atomselect top "all"]')
      df_tcl[6,1]<-paste0('$all set occupancy 0')
      df_tcl[7,1]<-paste0('set fixed [atomselect top " protein and resid  ',df_structure$C2[j],'"]')
      df_tcl[8,1]<-paste0('$fixed set occupancy 1')
      df_tcl[9,1]<-paste0('$all writepdb ',df_structure$bond[j],'_1.ref')
      df_tcl[10,1]<-paste0('mol delete all')
      write.table(df_tcl,paste0(v_start[i],'/protein/tcl/refgen_',df_structure$bond[j],'_1.tcl'),row.names = F,col.names = F,quote = F)
      df_tcl<-data.frame(matrix(nrow = 10,ncol = 1))
      df_tcl[1,1]<-paste0('cd ', part,v_start[i],'/protein')
      df_tcl[2,1]<-paste0('mol delete all')
      df_tcl[3,1]<-paste0('mol new ',df_structure$bond[j],'.psf')
      df_tcl[4,1]<-paste0('mol addfile ',df_structure$bond[j],'.pdb')
      df_tcl[5,1]<-paste0('set all [atomselect top "all"]')
      df_tcl[6,1]<-paste0('$all set occupancy 0')
      df_tcl[7,1]<-paste0('set fixed [atomselect top " protein and resid  ',df_structure$C1[j],'"]')
      df_tcl[8,1]<-paste0('$fixed set occupancy 1')
      df_tcl[9,1]<-paste0('$all writepdb ',df_structure$bond[j],'_2.ref')
      df_tcl[10,1]<-paste0('mol delete all')
      write.table(df_tcl,paste0(v_start[i],'/protein/tcl/refgen_',df_structure$bond[j],'_2.tcl'),row.names = F,col.names = F,quote = F) 
    }
  }
}
df_tcl<-data.frame(matrix(nrow=nrow(df_structure),ncol=length(v_start)))
for (i in 1:length(v_start)) {
  for (j in 1:nrow(df_structure)) {
    df_tcl[j,i]<-paste0("source ",part,v_start[i],'/protein/tcl/psfgen_',df_structure$bond[j],'.tcl\n',
                        "source ",part,v_start[i],'/protein/tcl/fixatom_',df_structure$bond[j],'_1.tcl\n',
                        "source ",part,v_start[i],'/protein/tcl/fixatom_',df_structure$bond[j],'_2.tcl\n',
                        "source ",part,v_start[i],'/protein/tcl/refgen_',df_structure$bond[j],'_1.tcl\n',
                        "source ",part,v_start[i],'/protein/tcl/refgen_',df_structure$bond[j],'_2.tcl\n')
    
  }
  df_tcl_out<-data.frame(matrix(nrow=1,ncol=length(v_start)))
  df_tcl_out[1,1]<-"exit now"
  colnames(df_tcl_out)<-colnames(df_tcl)
  df_tcl<-rbind(df_tcl,df_tcl_out)
  write.table(df_tcl,paste0(part,v_start[i],"/prepare_SMD_tcl.tcl"),row.names = F,col.names = F,quote = F,sep = "\n",na = "") 
  system(command = paste0("vmd -dispdev text -e ",part,v_start[i],"/prepare_SMD_tcl.tcl"),ignore.stdout=T,wait = T) 
}

