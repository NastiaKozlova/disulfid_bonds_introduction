part_start = commandArgs(trailingOnly=TRUE)

comand_to_run_namd<-"comand_to_run_namd "
library(bio3d)
library(ggplot2)
library(dplyr)
setwd(part_start)
number_y<-c(0:10000)
part<-paste0(part_start)
setwd(part)
df_structure<-read.csv("df_structure.csv",stringsAsFactors = F)
v_start<-c('SMD')

df_vstar<-data.frame(matrix(ncol=2,nrow = length(v_start)))
colnames(df_vstar)<-c("name","velosity")
df_vstar$name<-v_start
df_vstar$velosity<-c(0.01)
for (i in 1:nrow(df_vstar)) {
  l<-(1)
  if(!dir.exists(v_start[i])){dir.create(v_start[i])}
  if(!dir.exists(paste0(v_start[i],"/protein"))){dir.create(paste0(v_start[i],"/protein"))}
  if(!dir.exists(paste0(v_start[i],"/protein/tcl"))){dir.create(paste0(v_start[i],"/protein/tcl"))}
  if(!dir.exists(paste0(v_start[i],"/pdb"))){dir.create(paste0(v_start[i],"/pdb"))}
  if(!dir.exists(paste0(v_start[i],"/quench"))){dir.create(paste0(v_start[i],"/quench"))}
  if(!dir.exists(paste0(v_start[i],"/din"))){dir.create(paste0(v_start[i],"/din"))}
  for (j in 1:nrow(df_structure)) {
    if (file.exists(paste0("protein/start.pdb"))) {
      pdb<-read.pdb(paste0("protein/start.pdb"))
      df_pdb<-pdb$atom
      for (q in 1:2) {
        if (q==1){w<-2}else{w<-1}
        df_pdb_C<-df_pdb%>%filter(elety=="SG")  
        df_pdb_C1<-df_pdb_C%>%filter(resno==df_structure[j,(q+1)])
        df_pdb_C2<-df_pdb_C%>%filter(resno==df_structure[j,(w+1)])
        len<-sqrt((df_pdb_C1$x[1]-df_pdb_C2$x[1])^2+(df_pdb_C1$y[1]-df_pdb_C2$y[1])^2+(df_pdb_C1$z[1]-df_pdb_C2$z[1])^2)
        x<-(df_pdb_C1$x[1]-df_pdb_C2$x[1])/len
        y<-(df_pdb_C1$y[1]-df_pdb_C2$y[1])/len
        z<-(df_pdb_C1$z[1]-df_pdb_C2$z[1])/len
        
        df_tcl<-data.frame(matrix(nrow = 1,ncol = 1))
        df_tcl[1,1]<-paste0("#############################################################\n",
                            "## JOB DESCRIPTION                                         ##\n",
                            "#############################################################\n\n",
                            "# This is what this job does\n\n",
                            "#############################################################\n",
                            "## ADJUSTABLE PARAMETERS                                   ##\n",
                            "#############################################################\n\n",
                            "structure          protein/",df_structure$bonds[j],".psf\n",
                            "coordinates        protein/",df_structure$bonds[j],".pdb\n",
                            "#outputName         myoutput\n\n",
                            "set temperature    310\nfirsttimestep      0\n\n",
                            "#############################################################\n",
                            "## SIMULATION PARAMETERS                                   ##\n",
                            "#############################################################\n",
                            "# Input\n",
                            "paraTypeCharmm on\n",
                            "parameters ",part_start,"toppar/par_all36_prot.prm\n",
                            "parameters ",part_start,"toppar/par_all36_carb.prm\n",
                            "parameters ",part_start,"toppar/par_all36_cgenff.prm\n",
                            "parameters ",part_start,"toppar/par_all36_lipid.prm\n",
                            "parameters ",part_start,"toppar/par_all36_na.prm\n",
                            "parameters ",part_start,"toppar/toppar_water_ions_namd.str\n",
                            "# NOTE: Do not set the initial velocity temperature if you \n",
                            "# have also specified a .vel restart file!\n",
                            "temperature         $temperature\n",
                            "# Periodic Boundary conditions\n",
                            "# NOTE: Do not set the periodic cell basis if you have also \n",
                            "# specified an .xsc restart file!\n",
                            "if {0} { \n",
                            "cellBasisVector1    ",round(max(df_pdb$x)-min(df_pdb$x),digits = 0),"  0    0\n",
                            "cellBasisVector2     0   ",round(max(df_pdb$y)-min(df_pdb$y),digits = 0),"  0\n",
                            "cellBasisVector3     0    0   ",round(max(df_pdb$z)-min(df_pdb$z),digits = 0),"\n",
                            "cellOrigin           ",round(mean(df_pdb$x),digits = 0)," ", round(mean(df_pdb$y),digits = 0)," ", round(mean(df_pdb$z),digits = 0),"\n",
                            "}\n",
                            "wrapWater           on\n",
                            "wrapAll             on\n",
                            "# Force-Field Parameters\n",
                            "exclude             scaled1-4\n",
                            "1-4scaling          1.0\n",
                            "cutoff              12.0\n",
                            "switching           on\n",
                            "switchdist          10.0\n",
                            "pairlistdist        14.0\n",
                            "# Integrator Parameters\n",
                            "timestep            1.0  \n",
                            "rigidBonds          all  \n",
                            "nonbondedFreq       1\n",
                            "fullElectFrequency  2\n",
                            "stepspercycle       10\n",
                            #PME (for full-system periodic electrostatics)
                            #if {0} {
                            "#PME                 yes\n",
                            "#PMEGridSpacing      1.0\n",
                            
                            "#manual grid definition\n",
                            "PMEGridSizeX ",round(max(df_pdb$x)-min(df_pdb$x),digits = 0),"\n",
                            "PMEGridSizeY ",round(max(df_pdb$y)-min(df_pdb$y),digits = 0),"\n",
                            "PMEGridSizeZ ",round(max(df_pdb$z)-min(df_pdb$z),digits = 0),"\n",
                            #}
                            "restartfreq         500     ;# 500steps = every 1ps\n",
                            "dcdfreq             100\n",
                            "xstFreq             100\n",
                            "outputEnergies      100\n",
                            "outputPressure      100\n",
                            
                            "dcdfile quench/quench_",df_structure$bonds[j],"_",q,".dcd\n",
                            
                            "binaryoutput no\n",
                            "outputname pdb/",df_structure$bonds[j],"_",q,"\n",
                            
                            "#############################################################\n",
                            "## EXTRA PARAMETERS                                        ##\n",
                            "#############################################################\n",
                            
                            "# Put here any custom parameters that are specific to \n",
                            "# this job (e.g., SMD, TclForces, etc...)\n",
                            
                            "SMD	on\n",
                            "SMDFile protein/",df_structure$bonds[j],"_",q,".ref\n",
                            "SMDk	0.01\n",
 #                           "SMDVel	0.005\n",
                            "SMDVel	",df_vstar$velosity[i],"\n",
                            "SMDDir ",l*x," ",l*y," ",l*z,"\n",
                            "SMDOutputFreq 10000 \n",
                            
                            # Fixed Atoms Constraint (set PDB beta-column to 1)
                            "fixedAtoms          on\n",
                            "fixedAtomsFile      protein/",df_structure$bonds[j],"_",q,".fix\n",
                            "fixedAtomsCol       B\n",
                            "fixedAtomsForces    on\n",
                            "#############################################################\n",
                            "## EXECUTION SCRIPT                                        ##\n",
                            "#############################################################\n",
                            "# Minimization\n",
                            "if {0} {\n",
                            "minimize            100\n",
                            "reinitvels          $temperature\n",
                            "}\n",
                            "run 500000 ;# 0.5ns")
        write.table(df_tcl,paste0(v_start[i],"/",df_structure$bonds[j],"_",q,".conf"),col.names = F,row.names = F,quote = F)
      }
    }
  }
}
df_tcl<-data.frame(matrix(nrow=nrow(df_structure),ncol=1))
for (j in 1:nrow(df_structure)) {
  df_tcl[j,1]<-paste0("\n ",comand_to_run_namd," ",
                      "",df_structure$bonds[j],"_1.conf > ",df_structure$bonds[j],"_1.out",
                      "\n ",comand_to_run_namd," ",
                      "",df_structure$bonds[j],"_2.conf > ",df_structure$bonds[j],"_2.out")
}
df_tcl[1,1]<-paste0("cd ",part,v_start[1],"\n",df_tcl[1,1])
write.table(df_tcl,paste0(part,"namd_count.txt"),row.names = F,col.names = F,quote = F,sep = "\n")
