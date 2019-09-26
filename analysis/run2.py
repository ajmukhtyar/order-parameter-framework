LD_LIBRARY_PATH="/usr/local/intel/fc/10.0.023/lib:/fs/home/ajm529/lammps-11Aug17/src:/usr/local/mpich2/icse/lib"
import numpy as np
from lammps import lammps
lmp=lammps()
f1=open("checkcluster.txt","w",0)
lmp.file("control.dpd")
import os
import cluster3
import B_cluster
j=1000000
lmp.command("run 1000000")
lmp.command("write_data %s.data" %j)
os.system('./combinex.sh %s' %j)	
zz=cluster3.nmaxA(j)+B_cluster.nmaxB(j)
f1.write(str(zz)+"\n")

