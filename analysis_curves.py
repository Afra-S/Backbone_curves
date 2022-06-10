from glob import glob
import copy
import shutil
import subprocess
import os
import re
import string
import struct
import math
import numpy
import numpy as np
from scipy.stats import circmean
from scipy.stats import circstd
namefile=input("name file")

filetrj=open(namefile+'.pdb','r')
ki=-1
h,w=10000, 200;
parA=np.zeros((h, w,5))
parB=np.zeros((h, w,8))
parC=np.zeros((h, w,9))
puck=np.zeros((h, w))
pucktot=np.zeros((h, w,10))
puckhis=np.zeros((w,10))
dicpucker={"C1'en":0, "C1'ex":5,"C2'en":1,"C2'ex":6,"C3'en":2,"C3'ex":7,"C4'en":3,"C4'ex":8,"O1'en":4,"O1'ex":9}
# 0 0-1  9  1-2  17 2-3  24 3-4  30 4-5   35 5-6  39 6-7  42 7-8  44 8-9 
# 1 0-2  10 1-3  18 2-4  25 3-5  31 4-6   36 5-7  40 6-8  43 7-9
# 2 0-3  11 1-4  19 2-5  26 3-6  32 4-7   37 5-8  41 6-9
# 3 0-4  12 1-5  20 2-6  27 3-7  33 4-8   38 5-9
# 4 0-5  13 1-6  21 2-7  28 3-8  34 4-9
# 5 0-6  14 1-7  22 2-8  29 3-9
# 6 0-7  15 1-8  23 2-9 
# 7 0-8  16 1-9
# 8 0-9

pucktrans=np.zeros((w,45))

for line in filetrj:
    if(line[0:5]=='MODEL'):
        filesnap=open('snap.pdb','w')
        filesnap.write(line)
        ki=ki+1
        print(ki)
    elif(line[0:4]=='ATOM'):
        filesnap.write(line)
    elif(line[0:6]=='ENDMDL'):
        filesnap.write(line)
        filesnap.close()
        subprocess.call('sh curves_single.sh',shell=True)
        filec=open('bdna1.lis','r')
        kcon=-1
        log1=False
        log2=False
        log3=False
        ka=0
        for line2 in filec:
            check1=line2[0:11]
            check2=line2[2:13]
            check3=line2[2:14]
            check4=line2[2:25]
           # print(check1)
            if(check1=='  Strand  1'):
                nres=int(line2[15:20])
                #print(nres)
            if(log1==False):
                if(check2=='(A) BP-Axis'):
                    #print(line2)
                    log1=True
                    kl=0
            else:
                kl=kl+1
                if(kl>=2):
                    ka=ka+1
                    xdisp=float(line2[15:26])
                    ydisp=float(line2[27:35])
                    incl=float(line2[35:43])
                    tip=float(line2[43:53])
           
                    if(ka==1):
                        axbend=-1000
                    else:
                        axbend=float(line2[53:65])
                    
                    #print(line2[54:65])
                    #print(line2)
                    parA[ki,ka-1,0]=xdisp
                    parA[ki,ka-1,1]=ydisp
                    parA[ki,ka-1,2]=incl
                    parA[ki,ka-1,3]=tip
                    parA[ki,ka-1,4]=axbend
                    #print(parA[ki,ka-1,0:5])
                    #print(xdisp,ydisp,incl,tip,axbend,ka,nres)
                    #stop
                    if(ka==nres):
                        log1=False
            if(log2==False):
                kcon=kcon+1
               
                if(check3=='(C) Inter-BP'):
                    #print(line2)
                    log2=True
                    kl=0
                    ka=0
            else:
                kl=kl+1
                if(kl>=2):
                    ka=ka+1
                    #  print(line2[75:82])
                    #print(line2)
                    parB[ki,ka-1,0]=float(line2[18:26])
                    parB[ki,ka-1,1]=float(line2[26:36])
                    parB[ki,ka-1,2]=float(line2[36:44])
                    parB[ki,ka-1,3]=float(line2[44:52])
                    parB[ki,ka-1,4]=float(line2[51:59])
                    parB[ki,ka-1,5]=float(line2[59:67])
                    parB[ki,ka-1,6]=float(line2[67:75])
                    parB[ki,ka-1,7]=float(line2[75:82])
                  
                    
     
                    #if(kcon==280)
                    #print(parB[ki,ka-1,0:8])
                    if(ka==(nres-1)):
                        log2=False
                        
            if(log3==False):
                kcon=kcon+1
                #print(check4)
                if(check4=='(D) Backbone Parameters'):
                    print(line2)
                    log3=True
                    kl=0
                    ka=0
            else:
                kl=kl+1
                if(kl>=4):
                    ka=ka+1
                    #print(line2)
                    if(ka==1):
                        parC[ki,ka-1,0]=-1000
                        parC[ki,ka-1,1]=-1000
                    else:
                        parC[ki,ka-1,0]=float(line2[14:22])
                        parC[ki,ka-1,1]=float(line2[22:29])

                    pp=line2[79:84]
                    indxpuck=dicpucker[pp]
                    puck[ki,ka-1]=dicpucker[pp]
                    pucktot[ki,ka-1,indxpuck]=1
                    #print(puck[ki,ka-2],puck[ki,ka-1])
                    if(ki >0):
                        if(puck[ki-1,ka-1]!=puck[ki,ka-1]):
                            #print(ki,ka-1,puck[ki-1,ka-1],puck[ki,ka-1])
                             
                            if(puck[ki-1,ka-1] < puck[ki,ka-1]):
                                ixxx=int(puck[ki-1,ka-1])
                                iyyy= int(puck[ki,ka-1])
                            else:
                                iyyy=int(puck[ki-1,ka-1])
                                ixxx= int(puck[ki,ka-1])
                            if(ixxx==0):
                                indxmat=iyyy-1
                                pucktrans[ka-1,indxmat]=pucktrans[ka-1,indxmat]+1
                            else:
                                sumind=0
                                for ljj in range(1,ixxx+1):
                                    sumind=sumind+10-ljj
                                    #print(ixxx,iyyy,sumind)
                                indxmat=sumind+iyyy-ixxx-1
                                #print('fin',ixxx,iyyy,indxmat)
                                pucktrans[ka-1,indxmat]=pucktrans[ka-1,indxmat]+1
                    #print(dicpucker[pp])
                    parC[ki,ka-1,2]=float(line2[29:36])
                    parC[ki,ka-1,3]=float(line2[36:42])
                    puckhis[ka-1,dicpucker[pp]]=puckhis[ka-1,dicpucker[pp]]+1
                    parC[ki,ka-1,6]=float(line2[57:63])
                    parC[ki,ka-1,7]=float(line2[63:70])
                    parC[ki,ka-1,8]=float(line2[70:78])
                    if(ka==nres):
                        parC[ki,ka-1,4]=-1000
                        parC[ki,ka-1,5]=-1000
                        log3=False
                     #   stop
                    else:
                        parC[ki,ka-1,4]=float(line2[42:50])
                        parC[ki,ka-1,5]=float(line2[49:57])

kframes=ki+1
avera=np.zeros((nres,5))
fluca=np.zeros((nres,5))
hispuck=np.zeros((nres,10))
for i in range(0,nres):
    for j in range(0,2):
        avera[i,j]=np.mean(parA[0:kframes,i,j])
        fluca[i,j]=np.std(parA[0:kframes,i,j])
    for j in range(2,5):
        if(i==0 and j==4):
            avera[i,j]=-1000
            fluca[i,j]=-1000
        else: 
            avera[i,j]=circmean(parA[0:kframes,i,j]*math.pi/180,low=-math.pi,high=math.pi)*180/math.pi
            fluca[i,j]=circstd(parA[0:kframes,i,j]*math.pi/180,low=-math.pi,high=math.pi)*180/math.pi
        #print(i,avera[i,0],parA[0:kframes,i,0])
    #print(avera[i,2],parA[0:kframes,i,2])
    aaa=puckhis[i,0:10]/kframes
    for j in range(0,10):
        hispuck[i,j]=aaa[j]

np.savetxt('base_prop_ave.dat',avera,fmt='%10.2f')
np.savetxt('base_prop_fluc.dat',fluca,fmt='%10.2f')
np.savetxt('his_puck.dat',hispuck,fmt='%8.3f')

averb=np.zeros((nres-1,8))
flucb=np.zeros((nres-1,8))
for i in range(0,nres-1):
    for j in range(0,3):
        averb[i,j]=np.mean(parB[0:kframes,i,j])
        flucb[i,j]=np.std(parB[0:kframes,i,j])
    for j in range(3,8):
        averb[i,j]=circmean(parB[0:kframes,i,j]*math.pi/180,low=-math.pi,high=math.pi)*180/math.pi
        flucb[i,j]=circstd(parB[0:kframes,i,j]*math.pi/180,low=-math.pi,high=math.pi)*180/math.pi
   # print(averb[i,0:8])

np.savetxt('inter_base_ave.dat',averb,fmt='%10.2f')
np.savetxt('inter_base_fluc.dat',flucb,fmt='%10.2f')
    
averc=np.zeros((nres,9))
flucc=np.zeros((nres,9))
for i in range(0,nres):
    for j in range(0,8):
        if(i==0 and j<2):
            averc[i,j]=-1000
            flucc[i,j]=-1000
        elif(i==nres-1 and j==4):
            averc[i,j]=-1000
            flucc[i,j]=-1000
        elif(i==nres-1 and j==5):
            averc[i,j]=-1000
            flucc[i,j]=-1000
        else:
            averc[i,j]=circmean(parC[0:kframes,i,j]*math.pi/180,low=-math.pi,high=math.pi)*180/math.pi
            flucc[i,j]=circstd(parC[0:kframes,i,j]*math.pi/180,low=-math.pi,high=math.pi)*180/math.pi
        
    averc[i,8]=np.mean(parC[0:kframes,i,8])
    flucc[i,8]=np.std(parC[0:kframes,i,8])

np.savetxt('angular_conf_ave.dat',averc,fmt='%10.2f')
np.savetxt('angular_conf_fluc.dat',flucc,fmt='%10.2f')

matpuck=np.zeros((kframes,kframes))
for i in range(0,kframes):
    for j in range(i,kframes):
        dist=0
        for t in range(0,nres):
                if(puck[i,t]!=puck[j,t]):
                    dist=dist+1
        matpuck[i,j]=dist
        matpuck[j,i]=dist

        
np.savetxt('mat_puck.dat',matpuck,fmt='%4.2f')
np.savetxt('puck_trans.dat',pucktrans/(kframes-1),fmt='%5.3f')
