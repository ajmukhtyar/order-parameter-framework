#!/bin/bash


#LD_LIBRARY_PATH=/usr/local/intel/fc/10.0.023/lib:/usr/lib64:/usr/local/mpich2/icse/lib
LD_LIBRARY_PATH=/opt/ohpc/pub/software/vmd-1.9.3/lib64:/opt/ohpc/pub/mpi/openmpi3-gnu8/3.1.2/lib:/opt/ohpc/pub/compiler/gcc/8.2.0/lib64:/usr/lib64:/home/fs01/ajm529/lammps-16Mar18/src/

t0=$1

###Get Data####
sed -e s/SSSS/$t0/g getcoord.tcl > getcoordx.tcl
sed -e s/SSSS/$t0/g getcoordA.tcl > getcoordAx.tcl

module load vmd/1.9.3
vmd traj$t0.lammpstrj -dispdev text < getcoordx.tcl
vmd traj$t0.lammpstrj -dispdev text < getcoordAx.tcl

sed -e s/PPPP/$t0/g writedatafile.tcl > writedatafilex.tcl
vmd A$t0.lammpstrj -dispdev text < writedatafilex.tcl

sed -e s/PPPP/$t0/g B_writedatafile.tcl > B_writedatafilex.tcl
vmd B$t0.lammpstrj -dispdev text < B_writedatafilex.tcl

sed -e s/QQQQ/$t0/g controlq.dpd > controlx.dpd
/home/fs01/ajm529/lammps-16Mar18/src/lmp_serial -in controlx.dpd
sed -e s/QQQQ/$t0/g B_controlq.dpd > B_controlx.dpd
/home/fs01/ajm529/lammps-16Mar18/src/lmp_serial -in B_controlx.dpd


