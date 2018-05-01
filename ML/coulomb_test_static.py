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
MaxAtoms=10

MP_all=[]
MaterialsFile=pd.read_csv("test_static.csv")
for i in range(len(MaterialsFile)):
    MP_all.append(MaterialsFile['Materials ID'][i])
print(MP_all)
print(len(MP_all))

#test using small testing data
#Data=pd.read_csv('test.csv',index_col=False, sep=',')

f=open('CoulombResults_static.csv','w')
f.write('Material ID,')
for i in range(MaxAtoms):
    f.write('cm%s,' %str(i))
f.write('\n')

EValFiles=np.atleast_1d(['CoulombResults_static.csv'])
CMlists=[]

with MPRester(API_Key) as mp:# Materials Project API imported

    for i in range(0,len(MP_all)):
        #MPName=Data['Materials ID'][i]
        MPName=MP_all[i]
        StructDat = mp.query(criteria={"task_id": MPName}, properties=['structure','nsites'])
        #print(StructDat)
        nAtoms=int(StructDat[0]['nsites'])
        #print(nAtoms)
        print(MPName)
        #redefine the size of CMat
        MaxAtoms=nAtoms+1

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
    #WriteFiles(CMData['Materials ID'][i],CMData.iloc[i, -1],EVals,EValFiles)

f.close()
'''
print(CMdata)
df = pd.DataFrame(CMdata)
df.to_csv('CMdata.csv', index=False, header=False)'''
