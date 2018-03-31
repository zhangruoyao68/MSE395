'''
Created on Jan 6, 2018

@author: ethan
'''
import pandas as pd
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.core import periodic_table
import numpy as np
import pymatgen as mg
from scipy.linalg import eig
import csv

def WriteFiles(Chem,targets,descriptors,filelist):
    for iii in range(0,len(filelist)):
        f=open(filelist[iii],'a')
        f.write(Chem+',')
        for item in descriptors: #writes descriptors to file
            f.write(str(item)+',')
        f.write(str(targets)+'\n')
        f.close()
    return

API_Key='pF2RrzEPyZmIsLST'
MaxAtoms=200

CMData=pd.read_csv('test.csv')#input descriptor data file
#print(CMData["Materials ID"])
f=open('CoulombResults.csv','w')
f.write('Material ID,')
for i in range(0,MaxAtoms):
    f.write('cm%s,' % str(i))
f.write('\n')

EValFiles=np.atleast_1d(['CoulombResults.csv'])
CMlists=[]


with MPRester(API_Key) as mp:# Materials Project API imported

    for i in range(0,len(CMData)):
        MPName=CMData['Materials ID'][i]
        StructDat = mp.query(criteria={"task_id": MPName}, properties=['structure','nsites'])
        #print(StructDat)
        nAtoms=int(StructDat[0]['nsites'])
        #print(nAtoms)

        NuclearCharges=[]
        for qq in range(0,nAtoms):
            current=str(StructDat[0]['structure'].sites[qq]._species)
            elem=''.join(ii for ii in current if not ii.isdigit())
            NuclearCharges.append(mg.Element(elem).Z)

        CMat=np.zeros(MaxAtoms**2).reshape(MaxAtoms, MaxAtoms)
        #make offdiagonal components
        for k in range(0,nAtoms): #this makes lower triangle
            for kk in range(0,k):
                CMat[k][kk]=NuclearCharges[k]*NuclearCharges[kk]/(StructDat[0]['structure'].get_distance(k,kk))
        CMat=CMat+np.transpose(CMat)
        #make diagonal components
        for k in range(0,nAtoms):
            CMat[k][k]=0.5*NuclearCharges[k]**2.4

        #find eigenvals
        #print(CMat)
        vals,vecs=eig(CMat)
        #print(vals)
        #print(len(vals))
        EVals=[v.real for v in sorted(vals,reverse=True)]
        #CMlists.append(EVals)
        f.write(MPName)
        f.write(',')
        for item in EVals:
            f.write(str(item)+',')
        f.write('\n')
        print(MPName)
    #WriteFiles(CMData['Materials ID'][i],CMData.iloc[i, -1],EVals,EValFiles)

f.close()
'''
print(CMdata)
df = pd.DataFrame(CMdata)
df.to_csv('CMdata.csv', index=False, header=False)'''
