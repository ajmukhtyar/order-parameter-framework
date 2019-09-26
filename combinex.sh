#!/bin/bash

##Code to calculate the steinhardt bond order parameter Q6 using LAMMPS###
##Generates a dump file with the Q6 vector for every particle###

t0=$1
module load vmd/1.9.3


###For network A####
sed -e s/SSSS/$t0/g getcoordA.tcl > getcoordAx.tcl
vmd traj$t0.lammpstrj -dispdev text < getcoordAx.tcl
sed -e s/PPPP/$t0/g writedatafile.tcl > writedatafilex.tcl
vmd A$t0.lammpstrj -dispdev text < writedatafilex.tcl
sed -e s/QQQQ/$t0/g controlq.dpd > controlx.dpd
/home/fs01/ajm529/lammps-16Mar18/src/lmp_serial -in controlx.dpd

###For network B#####
sed -e s/SSSS/$t0/g getcoord.tcl > getcoordx.tcl
vmd traj$t0.lammpstrj -dispdev text < getcoordx.tcl
sed -e s/PPPP/$t0/g B_writedatafile.tcl > B_writedatafilex.tcl
vmd B$t0.lammpstrj -dispdev text < B_writedatafilex.tcl
sed -e s/QQQQ/$t0/g B_controlq.dpd > B_controlx.dpd
/home/fs01/ajm529/lammps-16Mar18/src/lmp_serial -in B_controlx.dpd


