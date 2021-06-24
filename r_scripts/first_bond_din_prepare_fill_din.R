#!/usr/bin/env R
part_start = commandArgs(trailingOnly=TRUE)
namd<-"comand_to_run_namd "
library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
number_y<-c(0:10000)
part<-paste0(part_start,'first_bond/')
setwd(part)
df_structure<-read.csv("min_length.csv",stringsAsFactors = F)

v_start<-c('SMD','SMD_rev')
if (!dir.exists('stablilisation')) {dir.create('stablilisation')}
if (!dir.exists('stablilisation/pdb')) {dir.create('stablilisation/pdb')}
if (!dir.exists('stablilisation/quench')) {dir.create('stablilisation/quench')}
if (!dir.exists('stablilisation/protein')) {dir.create('stablilisation/protein')}
i<-1
leng<-15
df_structure<-df_structure%>%mutate(scrpt=NA)
df_structure<-df_structure%>%mutate(namd=NA)
for (i in 1:nrow(df_structure)) {
  if (file.exists(paste0('protein_SMD/',df_structure$bonds[i],'.pdb'))) {
    pdb<-read.pdb(paste0('protein_SMD/',df_structure$bonds[i],'.pdb'))
    df_pdb<-pdb$atom
    write.pdb(pdb,paste0('stablilisation/protein/',df_structure$bonds[i],'.pdb'))
    fixed_1<-(df_structure$C1[i]-leng):(df_structure$C1[i]+leng)
    fixed_2<-(df_structure$C2[i]-leng):(df_structure$C2[i]+leng)
    if(min(fixed_2)-max(fixed_1)<leng){
      fixed_1<-(df_structure$C1[i]-leng):(df_structure$C2[i]-leng)
      fixed_2<-(df_structure$C1[i]+leng):(df_structure$C2[i]+leng)
    }
    fixed<-unique(c(fixed_1,fixed_2))
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0('#############################################################\n',
                        '## JOB DESCRIPTION                                         ##\n',
                        '#############################################################\n\n',
                        '# This is what this job does\n\n',
                        '#############################################################\n',
                        '## ADJUSTABLE PARAMETERS                                   ##\n',
                        '#############################################################\n\n',
                        'structure          protein/',df_structure$bonds[i],'.psf\n',
                        'coordinates        protein/',df_structure$bonds[i],'.pdb\n',
                        '#outputName         myoutput\n\n',
                        'set temperature    310\nfirsttimestep      0\n\n',
                        '#############################################################\n',
                        '## SIMULATION PARAMETERS                                   ##\n',
                        '#############################################################\n',
                        '# Input\n',
                        'paraTypeCharmm on\n',
                        'parameters ',part_start,'toppar/par_all36_prot.prm\n',
                        'parameters ',part_start,'toppar/par_all36_carb.prm\n',
                        'parameters ',part_start,'toppar/par_all36_cgenff.prm\n',
                        'parameters ',part_start,'toppar/par_all36_lipid.prm\n',
                        'parameters ',part_start,'toppar/par_all36_na.prm\n',
                        'parameters ',part_start,'toppar/toppar_water_ions_namd.str\n',
                        '# NOTE: Do not set the initial velocity temperature if you \n',
                        '# have also specified a .vel restart file!\n',
                        'temperature         $temperature\n',
                        '# Periodic Boundary conditions\n',
                        '# NOTE: Do not set the periodic cell basis if you have also \n',
                        '# specified an .xsc restart file!\n',
                        'if {0} { \n',
                        'cellBasisVector1    ',round(max(df_pdb$x)-min(df_pdb$x),digits = 0),'  0    0\n',
                        'cellBasisVector2     0   ',round(max(df_pdb$y)-min(df_pdb$y),digits = 0),'  0\n',
                        'cellBasisVector3     0    0   ',round(max(df_pdb$z)-min(df_pdb$z),digits = 0),'\n',
                        'cellOrigin           ',round(mean(df_pdb$x),digits = 0),' ', round(mean(df_pdb$y),digits = 0),' ', round(mean(df_pdb$z),digits = 0),'\n',
                        '}\n',
                        'wrapWater           on\n',
                        'wrapAll             on\n',
                        '# Force-Field Parameters\n',
                        'exclude             scaled1-4\n',
                        '1-4scaling          1.0\n',
                        'cutoff              12.0\n',
                        'switching           on\n',
                        'switchdist          10.0\n',
                        'pairlistdist        14.0\n',
                        '# Integrator Parameters\n',
                        'timestep            1.0  ;# 2fs/step\n',
                        'rigidBonds          all  ;# needed for 2fs steps\n',
                        'nonbondedFreq       1\n',
                        'fullElectFrequency  2\n',
                        'stepspercycle       10\n',
                        #PME (for full-system periodic electrostatics)
                        #if {0} {
                        '#PME                 yes\n',
                        '#PMEGridSpacing      1.0\n',
                        '#manual grid definition\n',
                        'PMEGridSizeX ',round(max(df_pdb$x)-min(df_pdb$x),digits = 0),'\n',
                        'PMEGridSizeY ',round(max(df_pdb$y)-min(df_pdb$y),digits = 0),'\n',
                        'PMEGridSizeZ ',round(max(df_pdb$z)-min(df_pdb$z),digits = 0),'\n',
                        #}
                        'restartfreq         500     ;# 500steps = every 1ps\n',
                        'dcdfreq             1000\n',
                        'xstFreq             1000\n',
                        'outputEnergies      1000\n',
                        'outputPressure      1000\n',
                        
                        'dcdfile quench/quench_',df_structure$bonds[i],'.dcd\n',
                        
                        'binaryoutput no\n',
                        'outputname pdb/',df_structure$bonds[i],'\n',
                        
                        '#############################################################\n',
                        '## EXTRA PARAMETERS                                        ##\n',
                        '#############################################################\n',
                        'fixedAtoms          on\n',
                        'fixedAtomsFile      protein/',df_structure$bonds[i],'.fix\n',
                        'fixedAtomsCol       B\n',
                        'fixedAtomsForces    on\n',
                        '#############################################################\n',
                        '## EXECUTION SCRIPT                                        ##\n',
                        '#############################################################\n',
                        '# Minimization\n',
                        'if {0} {\n',
                        'minimize            100\n',
                        'reinitvels          $temperature\n',
                        '}\n',
                        'run 500000 ;# 0.1ns')
    write.table(df_tcl,paste0('stablilisation/',df_structure$bonds[i],'.conf'),col.names = F,row.names = F,quote = F)
    df_structure$namd[i]<-paste0("cd ",part,"stablilisation/\n",
                                 namd, "",df_structure$bonds[i],".conf >","",df_structure$bonds[i],".out" )
    df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
    df_tcl[1,1]<-paste0('cd ',part,'stablilisation/','\n',
                        'mol delete all\n',
                        'package require psfgen\n', 
                        'resetpsf\n',
                        'topology ',part_start,'toppar/top_all36_prot.rtf\n',
                        'topology ',part_start,'toppar/toppar_water_ions_namd.str\n',
                        'pdbalias residue HIS HSE\n',
                        'pdbalias atom ILE CD1 CD\n',
                        'segment U { pdb protein/',df_structure$bonds[i],'.pdb\n','}\n',
                        'patch DISU U:',df_structure$C1[i],' U:',df_structure$C2[i],'\n',
                        'coordpdb protein/',df_structure$bonds[i],'.pdb U\n',
                        'regenerate angles dihedrals\n','guesscoord\n',
                        'guesscoord\n',
                        'writepdb protein/',df_structure$bonds[i],'.pdb\n',
                        'writepsf protein/',df_structure$bonds[i],'.psf\n',
                        'mol delete all\n',
                        'mol new protein/',df_structure$bonds[i],'.psf\n',
                        'mol addfile protein/',df_structure$bonds[i],'.pdb\n',
                        'set all [atomselect top "all"]\n',
                        '$all set beta 1\n',
                        'set fixed [atomselect top " protein and resid ',paste0(fixed,collapse = " "), '"]\n',
                        '$fixed set beta 0\n',
                        '$all writepdb protein/',df_structure$bonds[i],'.fix\n',
                        "\n\nexit now")
    write.table(df_tcl,paste0('stablilisation/protein/psfgen_',df_structure$bonds[i],'.tcl'),col.names = F,row.names = F,quote = F)
    system(command = paste0("vmd -dispdev text -e ", part,'stablilisation/protein/psfgen_',df_structure$bonds[i],'.tcl'),ignore.stdout=T,wait = T) 
  }
}
df_structure<-df_structure%>%select(namd)
df_structure<-df_structure%>%filter(!is.na(namd))
write.table(df_structure,paste0(part,"namd_stabilisation.txt"),col.names = F,row.names = F,quote = F)
