# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:54:25 2017

@author: ankitamukhtyar
"""

import sys
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

Natoms=13824
Lbox=23
frame=''
f1=open("coord%s.xyz" %frame)
index = np.zeros(Natoms)
x1=np.zeros(Natoms)
y1=np.zeros(Natoms)
z1=np.zeros(Natoms)
atomtype=np.zeros(Natoms)
avgdotprodts=np.zeros(Natoms)
sl=np.zeros([int(0.225*Natoms),3])
intcpt=np.zeros([Natoms,3])
rsq=np.zeros(Natoms)
clus=np.zeros(Natoms)
cyllike=np.zeros(Natoms)
mat=np.zeros([Natoms,Natoms])

##Reading coordinates###

p=0
for line in f1: 
    index[p]=float(line.split()[0])
    x=(float(line.split()[1]))
    y=(float(line.split()[2]))
    z=(float(line.split()[3]))
    if x<((-Lbox)*0.5): x = x+Lbox
    if x>=(Lbox*0.5): x = x-Lbox
    if y<((-Lbox)*0.5): y = y+Lbox
    if y>=(Lbox*0.5): y = y-Lbox
    if z<((-Lbox)*0.5): z = z+Lbox
    if z>=(Lbox*0.5): z = z-Lbox  
    x1[p]=x
    y1[p]=y
    z1[p]=z
    p=p+1
    
print "Got xyz"

###Finding cylinder signature vector for every particle by using Singular Value decomposition###

	
f3=open("sameneigh%s.txt" %frame)

def lineFit(points):
    datamean = points.mean(axis=0)
    uu, dd, vv = np.linalg.svd(points - datamean)
    dirvec=vv[0]
    return dirvec

f5=open("atomtype%s.txt" %frame)
p5=0
for line4 in f5:
	atomtype[p5]=float(line4.split()[1])
	p5=p5+1
 
 ### Finding highly aligned neighboring particles by measuring the dot product between their signature vectors ###

p2=0
for line2 in f3:
    if atomtype[p2]==1:
        n=len(line2.split()[1:])
        if n>0:
            distvec=np.zeros([3,n])
            r=np.zeros(n)
            p3=0
            for i in range(1,n+1):
                num=int(line2.split()[i])
                dx=(x1[p2]-x1[num])
                dy=(y1[p2]-y1[num])
                dz=(z1[p2]-z1[num])
                if abs(dx)>Lbox*0.5:
                    if dx>0: x1[num]=x1[num]+Lbox 
                    else:x1[num]=x1[num]-Lbox
                if abs(dy)>Lbox*0.5:
                    if dy>0: y1[num]=y1[num]+Lbox 
                    else:y1[num]=y1[num]-Lbox
                if abs(dz)>Lbox*0.5:
                    if dz>0: z1[num]=z1[num]+Lbox 
                    else:z1[num]=z1[num]-Lbox
                distvec[0][p3]=x1[num]
                distvec[1][p3]=y1[num]
                distvec[2][p3]=z1[num]
                p3=p3+1
            pts = np.concatenate((distvec[0][:, np.newaxis], distvec[1][:, np.newaxis], distvec[2][:, np.newaxis]), axis=1)
            sl[p2][0]=lineFit(pts)[0]
            sl[p2][1]=lineFit(pts)[1]
            sl[p2][2]=lineFit(pts)[2]
        p2=p2+1

print "Got line information"

## Clustering Algorithm for the lamellar like particles ####

p4=0
f4=open("neigh%s.txt" %frame)
count=0 
countdiff=0
countcyl=np.zeros(Natoms)   
for line3 in f4:
    if atomtype[p4]==1:
        n1=len(line3.split()[1:])
        ncyl=0
        cylind=[]
	diffatomtype=[]
        ####Finding cylinder like particles######
        for p in range(1,n1+1):
            nump=int(line3.split()[p])
            if nump!=p4:
                if atomtype[nump]==1:
                    dotprod1=np.dot(sl[p4],sl[nump])
                    if abs(dotprod1)>0.95:
                        ncyl=ncyl+1
                        cylind.append(nump)
		else: diffatomtype.append(nump)
        cyllike[p4]=ncyl
        if ncyl!=0: 
            cylind.append(p4)
        #########################################
        if len(cylind)>2:
            count=count+1 
	    countdiff=countdiff+len(diffatomtype)
            countcyl[p4]=1
    p4=p4+1

p5=0
f5=open("neigh%s.txt" %frame)
for line4 in f5:
    if atomtype[p5]==1:
        if countcyl[p5]==1:
            n2=len(line4.split()[1:])
            clusloc=[]
            cylind=[]
            for i in range(1,n2+1):
                nump=int(line4.split()[i])
                if atomtype[nump]==1:
                    if countcyl[nump]==1:
                        clusloc.append(clus[nump])
                        cylind.append(nump)
                else:
                    clusloc.append(clus[nump])
                    cylind.append(nump)
            maxclusloc=max(clusloc)
            if maxclusloc==0: 
                maxclus=max(clus)
                for i in range(len(cylind)):
                    num1=cylind[i]
                    clus[num1]=maxclus+1
            else:
                uniclusloc=list(set(clusloc))
                for i in range(len(cylind)):
                    num1=cylind[i]
                    clus[num1]= maxclusloc
                for k in range(len(uniclusloc)):
                    if uniclusloc[k]!=0:
                        clus=[maxclusloc if l==uniclusloc[k] else l for l in clus]
    p5=p5+1


print "Got cluster information"
print "Total number of ordered particles"
print count+countdiff
print "Minority component"
print count       
  

fid=open('clus%s.txt' %(frame),"w")
for i in range(len(clus)):
    fid.write(str(clus[i])+'\n')
fid.close()

          
counts,bins=np.histogram(clus,bins=2000)
findmaxclus=bins[counts.argsort()]
if len(findmaxclus)>1:
    if findmaxclus[-1]==0:
        max_X=findmaxclus[-2]
        print "Four largest clusters"
        print np.round(max_X),np.round(findmaxclus[-3]),np.round(findmaxclus[-4]),np.round(findmaxclus[-5])
        print "No of particles in each"
        print counts[counts.argsort()[-2]],counts[counts.argsort()[-3]],counts[counts.argsort()[-4]],counts[counts.argsort()[-5]]
        
    else:
        max_X=findmaxclus[-1]
        print "Four largest clusters"
        print np.round(max_X),np.round(findmaxclus[-2]),np.round(findmaxclus[-3]),np.round(findmaxclus[-4])
        print "No of particles in each"
        print counts[counts.argsort()[-1]],counts[counts.argsort()[-2]],counts[counts.argsort()[-3]],counts[counts.argsort()[-4]]


