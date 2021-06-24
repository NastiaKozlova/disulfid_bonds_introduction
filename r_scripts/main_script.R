part_start<-"/home/nastia/projects/current/NaPi2b/disulfid_bonds_introduction/"
#save_start_file to paste0(part_start,"protein/start.pdb")
#prepare conf and tcl files for for introduction of first disulfid bond
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/start_structure.R ",part_start),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_psfgen_SMD.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_prepre_conf_SMD.R ",part_start,'first_bond/'),ignore.stdout=T,wait = T)

#run namd2 
print(paste0("kate ",part_start,"first_bond/namd_count.txt"))
#system(command = paste0(part_start,"first_bond/namd_count.txt"),ignore.stdout=T,wait = T) 
#vmd SMD analysis
#uni
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_prepare_tcl_sencond_str.R ",part_start,'first_bond/'),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_pre_lenght_difulf_bonds_SMD.R ",part_start,'first_bond/'),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_lenght_difulf_bonds_SMD.R ",part_start,'first_bond/'),ignore.stdout=T,wait = T)
#to test
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_din_prepare_fill_din.R ",part_start),ignore.stdout=T,wait = T)
#namd
print(paste0("kate ",part_start,"first_bond/namd_stabilisation.txt"))
#system(command = paste0(part_start,"first_bond/namd_stabilisation.txt"),ignore.stdout=T,wait = T) 
##second
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_to_second_bond.R ",part_start),ignore.stdout=T,wait = T)
#fill din
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_bond_psfgen_SMD.R ",part_start),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_bond_prepre_conf_SMD.R ",part_start),ignore.stdout=T,wait = T)
#namd
#system(command = paste0(part_start,"second_bond/namd_count.txt"),ignore.stdout=T,wait = T) 
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_prepare_tcl_sencond_str.R ",part_start,'second_bond/'),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_pre_lenght_difulf_bonds_SMD.R ",part_start,'second_bond/'),ignore.stdout=T,wait = T)
system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/first_bond_lenght_difulf_bonds_SMD.R ",part_start,'second_bond/'),ignore.stdout=T,wait = T)

system(command = paste0("Rscript --vanilla  ",part_start,"r_scripts/second_bond_din_prepare_fill_din.R ",part_start),ignore.stdout=T,wait = T)

