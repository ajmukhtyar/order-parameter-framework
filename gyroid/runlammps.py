LD_LIBRARY_PATH="/usr/local/intel/fc/10.0.023/lib:/fs/home/ajm529/lammps-11Aug17/src:/usr/local/mpich2/icse/lib"
import numpy as np
from lammps import lammps
lmp=lammps()
import os
import cluster3
import B_cluster
j=50
f1=open("totord.txt","w",0)
f2=open("thermooutput.txt","w",0)
lmp.file("control2.dpd")

while j<10000000000:
	lmp.command("run 50")
	pe=lmp.get_thermo("pe")
	etotal=lmp.get_thermo("etotal")
	temp=lmp.get_thermo("temp")
	vol=lmp.get_thermo("vol")
	lmp.command("write_data %s.data" %j)
	os.system('./combinex.sh %s' %j)	
	#print("Starting cluster calculation")
	zz=cluster3.nmaxA(j)+B_cluster.nmaxB(j)
	f1.write(str(zz)+"\n")
	f2.write(str(pe)+"\t\t"+str(etotal)+"\t\t"+str(temp)+"\t\t"+str(vol)+"\n")
	if zz<8000:
		os.system('rm atomtype%s.txt sameneigh%s.txt dump%s.xyz A%s.lammpstrj A%s.data' %(j,j,j,j,j))
		os.system('rm B_atomtype%s.txt B_sameneigh%s.txt B_dump%s.xyz B%s.lammpstrj B%s.data' %(j,j,j,j,j))
		j=j+50
		continue
		
	else:
		break

