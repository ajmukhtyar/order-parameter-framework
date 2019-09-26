# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 17:54:25 2017

@author: ankitamukhtyar
"""

def nmaxB(frame):
    import numpy as np
    Natoms=6912
    f1=open("B_dump%s.xyz" %frame)
    index = np.zeros(Natoms)
    Qlm_real=[[] for i in range(Natoms)]
    Qlm_im=[[] for i in range(Natoms)]
    Q6=np.zeros(Natoms)
    
    for line in f1:
        if line.split()[-1]=='c_myBond[27]':
            break
    p=0
    for line in f1: 
        index[p]=float(line.split()[0])
        Q6[p]=(float(line.split()[2]))
        n=len(line.split()[3:])
        for i in range(n):
            num=(float(line.split()[3+i]))
            if i%2==0:
                Qlm_real[p].append(num)
            else:
                Qlm_im[p].append(num)
        p=p+1
        
    atomtype=np.zeros(Natoms)
    f2=open("B_atomtype%s.txt" %frame)
    p1=0
    for line1 in f2:
        atomtype[p1]=float(line1.split()[1])
        p1=p1+1
    
    f3=open("B_sameneigh%s.txt" %frame)
    p2=0
    nlist=[[] for i in range(Natoms)]
    for line2 in f3: 
        if atomtype[p2]==2:
            n=len(line2.split()[1:])
            for j in range(1,n+1):
                num = float(line2.split()[j])
                nlist[p2].append(num)
            p2=p2+1
    
    dotprod=[[]for i in range(Natoms)]
    
    def Qlm(self,neigh):
        sum0=0
        for j in range(len(Qlm_real[self])):
            dr=Qlm_real[self][j]*Qlm_real[neigh][j]
            di=Qlm_im[self][j]*Qlm_im[neigh][j]
            drdi=dr+di
            sum1=sum0+drdi
            sum0=sum1
        return sum0
            
    def dotprod(self,neigh):
         numer=Qlm(self,neigh)
         denom=np.sqrt(Qlm(self,self)*Qlm(neigh,neigh))
         d=numer/denom
         return d
         
    Natoms=6912
    atomtype=np.zeros(Natoms)
    clus=np.zeros(int(Natoms))
    gyrlike=np.zeros(int(Natoms))
    output=open("B_output%s.txt" %frame,"w",0)
    f1=open("B_atomtype%s.txt" %frame)
    p1=0
    for line1 in f1:
        atomtype[p1]=float(line1.split()[1])
        p1=p1+1
    p4=0
    f2=open("B_sameneigh%s.txt" %frame)
    countgyr=np.zeros(Natoms)
    count=0
    for line3 in f2:
        if atomtype[p4]==2:
            n1=len(line3.split()[1:])
            clusloc=np.zeros(n1)
            ngyr=0
            gyrind=[]
            ####Finding gyroid like particles######
            for p2 in range(1,n1+1):
                nump=int(line3.split()[p2])
                if nump>=6912:nump=nump-6912
                dp=dotprod(p4,nump)
                if abs(dp)>0.7:
                    ngyr=ngyr+1
                    gyrind.append(nump)
            gyrlike[p4]=ngyr
            if ngyr!=0: gyrind.append(p4)
            #if p4==70: gyrindsave=gyrind
            #########################################
            if gyrlike[p4]>5:
                countgyr[p4]=1
                count=count+1
            p4=p4+1
    
    f3=open("B_sameneigh%s.txt" %frame)
    p5=0
    #tryx=1
    for line4 in f3:
        if atomtype[p5]==2:
            if countgyr[p5]==1:
                lamind=[]
                n2=len(line4.split()[1:])
                clusloc=[]
                for i in range(1,n2+1):
                    nump=int(line4.split()[i])
                    if nump>=6912:nump=nump-6912
                    if countgyr[nump]==1:
                        clusloc.append(clus[nump])
                        lamind.append(nump)
                if clusloc!=[]:
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
    
    output.write("Got cluster information"+"\n")
    output.write("Total ordered particles"+"\n")
    output.write(str(count)+"\n")
    
    
    fid=open('B_clus%s.txt' %(frame),"w")
    for i in range(len(clus)):
        fid.write(str(clus[i])+'\n')
    fid.close()
    
              
    counts,bins=np.histogram(clus,bins=2000)
    findmaxclus=bins[counts.argsort()]
    if len(findmaxclus)>1:
        if findmaxclus[-1]==0:
            max_X=findmaxclus[-2]
            output.write("Three largest clusters"+"\n")
            output.write(str(np.round(max_X))+ " " +str(np.round(findmaxclus[-3]))+" "+str(np.round(findmaxclus[-4]))+"\n")
            output.write("No of particles in each"+"\n")
            output.write(str(counts[counts.argsort()[-2]])+" "+str(counts[counts.argsort()[-3]])+" "+str(counts[counts.argsort()[-4]])+"\n")
            size=counts[counts.argsort()[-2]]
        else:
            max_X=findmaxclus[-1]
            output.write("Three largest clusters"+"\n")
            output.write(str(np.round(max_X))+ " " +str(np.round(findmaxclus[-2]))+" "+str(np.round(findmaxclus[-3]))+"\n")
            output.write("No of particles in each"+"\n")
            output.write(str(counts[counts.argsort()[-1]])+" "+str(counts[counts.argsort()[-2]])+" "+str(counts[counts.argsort()[-3]])+"\n")
            size=counts[counts.argsort()[-1]]
    return size

if __name__ == "__main__":
        nmaxB()
