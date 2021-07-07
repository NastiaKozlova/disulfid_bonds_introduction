# disulfid_bonds_introduction
Scripts to introduce disulfide bonds in transmembrane protein structure

Folder content:
protein - put original protein stucture
toppar - put CHARMM Force fields
r_scripts -folder with scripts

r_scripts/main_script.R - main script from which all scrips are rinning
r_scripts/start_structure.R - find disulfide bonds and divide them into three groups by z coordinate of CA atom of cysteins
TMD=cysteins in membrane
cyto=with z<(-18)
extra= cysteins with z>(18)
Separation of cysteins in heavily dependent on z coordinates, therefore You should be wery careful with this script r_scripts/start_structure.R and check it wery carefully or mofiy it if needed.
You should run namd 4 times (scripts for namd are in
first_bond/namd_count.txt
first_bond/namd_stabilisation.txt
second_bond/namd_count.txt
second_bond/namd_stabilisation.txt)
 
r_scripts/main_script.R - main script from which all scrips are rinning
r_scripts/start_structure.R - find disulfide bonds and divide them into three groups by z coordinate of CA atom of cysteins TMD=cysteins in membrane cyto=with z<(-18) extra= cysteins with z>(18) Separation of cysteins in heavily dependent on z coordinates, therefore You should asiste carefully output file, from the midle of the file and adgaste accordingly:"protein/cysteins.csv"
 
r_scripts/start_structure.R and check it wery carefully or mofiy it if needed. You should run namd 4 times (scripts for namd are:
first_bond/namd_count.txt
first_bond/namd_stabilisation.txt
second_bond/namd_count.txt
second_bond/namd_stabilisation.txt)
 
*Scripts are work in 4 major parts:
+ preparation systems to SMD(Steered molecular dynamics) MD (R scripts)
+ SMD MD (NAMD)
+ Find all possible disulfide bonds and prepare systems to stabilisation MD (R scripts)
+ Stabilisation MD (NAMD)
+ find all stable disulfide bonds... and and start step 1 again
 
 
All MD simulations are in vacuum.
 
All aminoacids exsept 10 around cysteins wich forme disulfide bonds are fixed during SMD, if length between CYS less then 30 AA then all aminoacids between CYS are also flexible. All aminoacids exsept 15 around cysteins wich formed disulfide bonds are fixed, if length between CYS less then 45 AA then all aminoacids between CYS are also flexible. Force to create disulfide bond is from one S of the CYS pair to another and from second to first. Force with with we pull S atoms are 0.01 kcal/mol/Å$ ^2$
Length of SMD and stabilased MD are 0.5ns.
 
Disulfide bond is counting as stable if system is stable during stabilasing MD are 0.5ns.
 