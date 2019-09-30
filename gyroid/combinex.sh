##Calculate bond order parameter Q6 using LAMMPS##
##Returns a dump file with the Q6 vector for every particle##

#!/bin/bash
module load vmd

### Get data for network A ####
vmd traj.lammpstrj -dispdev text < getcoordA.tcl
vmd A.lammpstrj -dispdev text < A_writedatafile.tcl
/home/fs01/ajm529/lammps-16Mar18/src/lmp_serial -in A_controlq.dpd

### Get data for network B ####
vmd traj.lammpstrj -dispdev text < getcoordB.tcl
vmd B.lammpstrj -dispdev text < B_writedatafile.tcl
/home/fs01/ajm529/lammps-16Mar18/src/lmp_serial -in B_controlq.dpd


