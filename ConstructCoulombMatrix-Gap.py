'''
Created on Jan 6, 2018

@author: ethan
'''
import pandas as pd
from pymatgen.matproj.rest import MPRester, MPRestError
from pymatgen.core import periodic_table
import numpy as np
import pymatgen as mg
from scipy.linalg import eig

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
MaxAtoms=160

PlasmonicData=pd.read_csv('CompoundDescriptors-Gap.csv',index_col=False, sep=',')#input descriptor data file
#print(PlasmonicData)

EValFiles=np.atleast_1d(['CoulombEValsGap.csv'])

with MPRester(API_Key) as mp:# Materials Project API imported

    for i in range(0,len(PlasmonicData)):
        print(i)
        MPName='mp-'+str(PlasmonicData['MPID'][i])
        StructDat = mp.query(criteria={"task_id": MPName}, properties=['structure','nsites'])
        nAtoms=int(StructDat[0]['nsites'])
        
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
        vals,vecs=eig(CMat)
        EVals=[v.real for v in sorted(vals,reverse=True)]
        WriteFiles(PlasmonicData['Material'][i],PlasmonicData.iloc[i, -1],EVals,EValFiles)
        
        #find max row
        #rowvals=CMat.sum(axis=0)
        #maxrow=np.argmax(rowvals)
        #MaxRow=sorted(CMat[maxrow],reverse=True)
        #WriteFiles(PlasmonicData['Chemical'][i],PlasmonicData.iloc[i, -1].values.tolist(),MaxRow,RowFiles)