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