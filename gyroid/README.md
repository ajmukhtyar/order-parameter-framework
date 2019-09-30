# Gyroid Phase Order Parameter 

Code uses Python v2.7 and VMD to run

For any given LAMMPS trajectory, the number of ordered (gyroid) particles can be found by running the following commands

1. ./combinex.sh

This bash script uses VMD to return information about the coordinates, atom type and nearest neighbors of every particle in both networks A and B and stores this information in files - atomtypeA/B.txt and sameneighA/B.txt

It also calculate the Q6 bond order parameter for every particle through files A/B_controlq.dpd which gets stored in the respective dump files. 

2. python A/B_cluster.py 

This calculates the gyroid signature vector for every particle and uses a clustering algorithm to find the size of the largest ordered nucleus. Information is stored in output. 
