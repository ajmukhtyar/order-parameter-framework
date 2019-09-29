Code uses Python v2.7 and VMD to run

For any given LAMMPS trajectory, the number of ordered (lamellar) particles can be found by running the following commands

1. vmd traj.lammpstrj -dispdev text < getcoord.tcl 

This uses VMD to return information about the coordinates, atom type and nearest neighbors of every particle and stores this information in files - coord.xyz, atomtype.txt, neigh.txt, and sameneigh.txt

2. python planarOP.py > output

This calculates the lamellar signature vector for every particle and uses a clustering algorithm to find the size of the largest ordered nucleus. 
Information is stored in output. 

The clusters can be visualized in VMD by loading the trajectory and running "source readclus.tcl" in the VMD terminal. The output file contains the index of the largest cluster (in this case 16) which can be seen by typing "beta = 16" in the graphical representation box in VMD.