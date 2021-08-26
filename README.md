# disulfid_bonds_introduction
Scripts to introduce disulfide bonds to the transmembrane protein structures.

Folder content:

protein - put original protein structure

toppar - put CHARMM Force fields

r_scripts -folder with scripts

r_scripts/main_script.R - main script from which all scrips are rinning

r_scripts/start_structure.R - find disulfide bonds and divide them into three groups by z coordinate of CA atom of cysteines

TMD=cysteins in membrane
cyto=with z<(-18)
extra= cysteines with z>(18)

The separation of cysteines is heavily dependent on z coordinates. You should be very careful with this script r_scripts/start_structure.R and check it very carefully modify it if needed.

You should run namd four times. Scripts for namd are in

first_bond/namd_count.txt
first_bond/namd_stabilisation.txt
second_bond/namd_count.txt
second_bond/namd_stabilisation.txt)
 
r_scripts/main_script.R - main script from which all scrips are rinning
r_scripts/start_structure.R - find disulfide bonds and divide them into three groups by z coordinate of CA atom of cysteines TMD=cysteins in membrane cyto=with z<(-18) extra= cysteines with z>(18) Separation of cysteines heavily dependent on z coordinates, therefore You should verify output file, from the middle of the file and adjust accordingly: "protein/cysteins.csv"
 
*Scripts are work in 4 major parts:
+ preparation systems to SMD(Steered molecular dynamics) MD (R scripts)
+ SMD MD (NAMD)
+ Find all possible disulfide bonds and prepare systems to stabilisation MD (R scripts)
+ Stabilisation MD (NAMD)
+ find all stable disulfide bonds... and start step 1 again
 
 
All MD simulations are in a vacuum.
 
All amino acids except 10 around cysteines which form disulfide bonds are fixed during SMD if the length between CYS less than 30 AA then all aminoacids between CYS are also flexible. All amino acids except 15 around cysteines which formed disulfide bonds are fixed if the length between CYS less than 45 AA then all aminoacids between CYS are also flexible. Force to create a disulfide bond is from one S of the CYS pair to another and from second to first. Force with we pull S atoms are 0.01 kcal/mol/Å$ ^2$
Length of SMD and stabilised MD are 0.5ns.
 
The disulfide bond is counting as stable if the system is stable during stabilising MD are 0.5ns.