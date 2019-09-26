
import numpy as np
from lammps import lammps
lmp=lammps()
import os
import A_cluster
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
	zz=A_cluster.nmaxA(j)+B_cluster.nmaxB(j)
	f1.write(str(zz)+"\n")
	f2.write(str(pe)+"\t\t"+str(etotal)+"\t\t"+str(temp)+"\t\t"+str(vol)+"\n")
	if zz<8000:
		j=j+50
		continue
		
	else:
		break

