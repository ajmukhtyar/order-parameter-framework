# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:54:25 2017

@author: ankitamukhtyar
"""

import sys
import numpy as np

Natoms=13824
Lbox=22
frame=""
f1=open("coord%s.xyz" %frame)
index = np.zeros(Natoms)
x1=np.zeros(Natoms)
y1=np.zeros(Natoms)
z1=np.zeros(Natoms)
atomtype=np.zeros(Natoms)
normvec=np.zeros([Natoms,3])
clus=np.zeros(Natoms)
lamlike=np.zeros(Natoms)


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


###Finding lamellar signature vector for every particle by using Singular Value decomposition###

f3=open("sameneigh%s.txt" %frame)

def planeFit(points):
    datamean = points.mean(axis=0)
    uu, dd, vv = np.linalg.svd(points - datamean)
    dirvec=vv[-1]
    return dirvec

p2=0
for line2 in f3:
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
        pts = np.concatenate((distvec[0][:, np.newaxis], distvec[1][:, np.newaxis], distvec[2][:, np.newaxis]),axis=1)
        normal = [planeFit(pts)[0],planeFit(pts)[1],planeFit(pts)[2]]
        normvec[p2][0]=normal[0]
        normvec[p2][1]=normal[1]
        normvec[p2][2]=normal[2]
    p2=p2+1
        

print "Got normal vectors"


### Finding highly aligned neighboring particles by measuring the dot product between their signature vectors ###
p4=0
f4=open("neigh%s.txt" %frame)
t=0
count=0
countlam=np.zeros(Natoms)
for line3 in f4:
    n1=len(line3.split()[1:])
    clusloc=np.zeros(n1)
    nlam=0
    lamind=[]
    ####Finding lamellar like particles######
    for p in range(1,n1+1):
        nump=int(line3.split()[p])
        if nump!=p4:
            dotprod1=np.dot(normvec[p4],normvec[nump])
            if abs(dotprod1)>0.95:
                nlam=nlam+1
                lamind.append(nump)
    lamlike[p4]=nlam
    if nlam!=0: lamind.append(p4)
    #########################################
    if lamlike[p4]>5:
        count=count+1
        countlam[p4]=1
    p4=p4+1


## Clustering Algorithm for the lamellar like particles ####
p5=0
f5=open("neigh%s.txt" %frame)
for line4 in f5:
    if countlam[p5]==1:
        lamind=[]
        n2=len(line4.split()[1:])
        clusloc=[]
        for i in range(1,n2+1):
            nump=int(line4.split()[i])
            if countlam[nump]==1:
                clusloc.append(clus[nump])
                lamind.append(nump)
        maxclusloc=max(clusloc)
        if maxclusloc==0: 
            maxclus=max(clus)
            for j in range(len(lamind)):
                num1=lamind[j]
                clus[num1]=maxclus+1
        else:
            uniclusloc=list(set(clusloc))
            for j in range(len(lamind)):
                num1=lamind[j]
                clus[num1]= maxclusloc
            for k in range(len(uniclusloc)):
                if uniclusloc[k]!=0:
                    clus=[maxclusloc if l==uniclusloc[k] else l for l in clus]
    p5=p5+1

print "Got cluster information"
print "Total ordered particles"
print count

fid=open('clus%s.txt' %frame,"w")
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
        print counts[counts.argsort()[-1]],counts[counts.argsort()[-2]],counts[counts.argsort()[-3]], counts[counts.argsort()[-4]] 

