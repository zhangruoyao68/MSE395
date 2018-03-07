'''
Created on May 14, 2017

@author: ethan
'''
#Makes descriptor file
#Includes DFT bandgap from Materials Project
import random
import os, glob
import shutil
import itertools
import sys
import re
import numpy as np
import math
import fileinput
import mendeleev
import random
import subprocess
import pandas as pd
from mendeleev import element
from numpy import genfromtxt
from tempfile import mkstemp
from shutil import move
from os import remove, close
from itertools import combinations
import numpy as np
import time
from pymatgen.ext.matproj import MPRester, MPRestError

start_time = time.time()

def SymmGroupFinder(GroupNum):
    if 1<=GroupNum<=2:
        Group='Triclinic'
    elif 3<=GroupNum<=15:
        Group='Monoclinic'
    elif 16<=GroupNum<=74:
        Group='Orthorhombic'
    elif 75<=GroupNum<=142:
        Group='Tetragonal'
    elif 143<=GroupNum<=167:
        Group='Trigonal'
    elif 168<=GroupNum<=194:
        Group='Hexagonal'
    elif 195<=GroupNum<=230:
        Group='Cubic'
    else:
        Group='Amorphous'
    return Group


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def NameAnalyzer(NameToTest):
    NameToTest= [ y for y in list(itertools.chain(*[re.split(r'\"(.*)\"', x) 
        for x in re.split(r'\((.*)\)', NameToTest)])) 
        if y != ''] #this splits at parentheses
    
    ParenSeparatedName=[]
    for j in range(0,len(NameToTest)):#This will check if first character is a number, due to parentheses
        firstChar=NameToTest[j][0]
        if firstChar.isdigit()==True:
            ParenSeparatedName.extend([a for a in re.split(r'([A-Z][a-z]*\d*)', NameToTest[j]) if a])
        else:
            #print NameToTest[j]
            ParenSeparatedName.append(NameToTest[j])
    
    SeparatedName=[]
    for j in range(0,len(ParenSeparatedName)):
        TempSegment=[a for a in re.split(r'([A-Z][a-z]*)', ParenSeparatedName[j]) if a]
        multiplier=1.0
        if isfloat(TempSegment[0])==False:
            if j<len(ParenSeparatedName)-1 and isfloat(ParenSeparatedName[j+1])==True:
                multiplier=float(ParenSeparatedName[j+1])
            for k in range(0,len(TempSegment)):
                #print(TempSegment[k],TempSegment[k].isnumeric())
                if isfloat(TempSegment[k])==False:
                    if k<len(TempSegment)-1 and isfloat(TempSegment[k+1])==True:
                        SeparatedName.append(TempSegment[k])
                        SeparatedName.append(str(float(TempSegment[k+1])*multiplier))
                    elif k<len(TempSegment)-1 and isfloat(TempSegment[k+1])==False:
                        SeparatedName.append(TempSegment[k])
                        SeparatedName.append(str(multiplier))
                    elif k==len(TempSegment)-1 and isfloat(TempSegment[k])==False:
                        SeparatedName.append(TempSegment[k])
                        SeparatedName.append(str(multiplier))
    return(SeparatedName)





#MaterialsFile=pd.read_csv("MaterialBandGaps.csv")
#MaterialsFile=pd.read_csv("test.csv")
MaterialsFile=np.genfromtxt('CNData.txt')
#print(MaterialsFile)
#print(MaterialsFile[0])
f=open('CompoundDescriptorsMP.csv','w')
#print(MaterialsFile)
Compounds=[]
#TcValues=[]
for i in range(0,len(MaterialsFile)):
    print(i)
    Compounds.append(int(MaterialsFile[i][0]))
    #TcValues.append(MaterialsFile[i][1])
print(Compounds)
numMaterials=len(Compounds)

#Set up datafile to write
#Eneg is maximum difference
#f.write('Material,sTot,pTot,dTot,fTot,RedMass,AverageMass,DevMass,AvgENeg,DiffENeg,DevENeg,OmegaAvg,OmegaDev,OmegaMax,OmegaMin,AvgIonicRad,DevIonicRad,AvgDipole,DevDipole,DiffDipole,AvgBlock,DiffBlock,TotalVolume,totalAtoms'+'\n')

ElementList=np.genfromtxt('NameListFile.txt',dtype=None)
Elements=[]
for i in range(0,len(ElementList)):
    Elements.append(ElementList[i][1].decode("utf-8"))
ElemData=[]#[[]]*len(Elements)
block = {'s':0,'p':1,'d':2,'f':3}
for elem in Elements:
    row=[elem]
    sTot=0
    pTot=0
    dTot=0
    fTot=0
    elconfig=(element(elem).econf)
    sEl=re.findall(r'(?<=s)(\d*)',elconfig)
    pEl=re.findall(r'(?<=p)(\d*)',elconfig)
    dEl=re.findall(r'(?<=d)(\d*)',elconfig)
    fEl=re.findall(r'(?<=f)(\d*)',elconfig)
    #print(elconfig,sEl,pEl,dEl,fEl)
    if sEl!=[]:
        sTot+=int(sEl[0])
    if pEl!=[]:
        pTot+=int(pEl[0])
    if dEl!=[]:
        dTot+=int(dEl[0])
    if fEl!=[]:
        fTot+=int(fEl[0])
    row.extend([sTot,pTot,dTot,fTot])
    row.append(element(elem).mass)
    row.append(element(elem).en_pauling)
    row.append(element(elem).ionenergies.get(1))
    row.append(element(elem).atomic_radius)
    row.append(element(elem).ionic_radii[0].ionic_radius)
    row.append(element(elem).dipole_polarizability)
    row.append(float(block[element(elem).block]))
    row.append(element(elem).atomic_radius**3*10**-7)
    ElemData.append(row)
    #print(row)
#print(ElemData)

ElementData={d[0]: d[1:] for d in ElemData}

DOSData=np.genfromtxt('IntegratedDos.csv',delimiter=',',dtype=None)
for i in range(0,len(DOSData)):
    Element=DOSData[i][0].decode("utf-8").split("_")[0]
    if Element in ElementData:
        #print(Element,ElementData[Element])
        ElementData[Element].append(10**5*DOSData[i][4])
    
suscData=np.genfromtxt('MagneticSusceptibility.csv',delimiter=',',dtype=None)
for i in range(0,len(suscData)):
    Element=suscData[i][0].decode("utf-8").split("_")[0]
    if Element in ElementData:
        ElementData[Element].append(suscData[i][1])

API_Key='pF2RrzEPyZmIsLST'
m=MPRester(API_Key)


for i in range(0,numMaterials):
    #Comp,MPID=Compounds[i].split("_")
    #Comp=Compounds[i]
    #print(Comp)
    #MatName=NameAnalyzer(Comp)
    #if 'Ac' in MatName:
    #    continue
    stabdata = m.query(criteria={'material_id':'mp-'+str(Compounds[i])}, properties=["pretty_formula", "e_above_hull",
                                                                      "formation_energy_per_atom","band_gap",
                                                                      "spacegroup.number",'volume',
                                                                      'structure','density','e_above_hull','nsites'])
    
    #print(stabdata)
    #print(stabdata[0]['spacegroup.number'])
    #stabdata = m.query(criteria={"pretty_formula": Comp}, properties=["material_id"])
    #if stabdata:
        #mini=min(stabdata, key=lambda x:x['formation_energy_per_atom'])
        #if 'material_id' in mini:
    symm=SymmGroupFinder(stabdata[0]['spacegroup.number'])
    num = Compounds[i]
    gap = stabdata[0]['band_gap']
    volume=str.format("{0:.6f}",stabdata[0]['volume'])
    formE=str.format("{0:.6f}",stabdata[0]['formation_energy_per_atom'])
    lattA=str.format("{0:.6f}",stabdata[0]['structure'].lattice.a)
    lattB=str.format("{0:.6f}",stabdata[0]['structure'].lattice.b)
    lattC=str.format("{0:.6f}",stabdata[0]['structure'].lattice.c)
    lattAlpha=str.format("{0:.6f}",stabdata[0]['structure'].lattice.angles[0])
    lattBeta=str.format("{0:.6f}",stabdata[0]['structure'].lattice.angles[1])
    lattGamma=str.format("{0:.6f}",stabdata[0]['structure'].lattice.angles[2])
    density=str.format("{0:.6f}",stabdata[0]['density'])
    eHull=str.format("{0:.6f}",stabdata[0]['e_above_hull'])
    NAtom=str(int(stabdata[0]['nsites']))
    Comp=stabdata[0]['pretty_formula']
    MatName=NameAnalyzer(Comp)
    #if ('Ac' or 'Pu') in MatName:
    #    continue
    print(Comp,num,gap)
    
    try:
        AtomicWeights=[float(k) for k in MatName[1::2]]
        AtomicSpecies=[k for k in MatName[0::2]]
        totalAtoms=sum([float(k) for k in MatName[1::2]]) #total number of atoms in formula unit
        #print(MatName,AtomicSpecies,AtomicWeights,totalAtoms)
        sTot=0
        pTot=0
        dTot=0
        fTot=0
        for j in range(0,len(AtomicSpecies)): #this gets reduced electronic configuration
            sTot+=ElementData[AtomicSpecies[j]][0]*AtomicWeights[j]
            pTot+=ElementData[AtomicSpecies[j]][1]*AtomicWeights[j]
            dTot+=ElementData[AtomicSpecies[j]][2]*AtomicWeights[j]
            fTot+=ElementData[AtomicSpecies[j]][3]*AtomicWeights[j]
        sTot=sTot/totalAtoms
        pTot=pTot/totalAtoms
        dTot=dTot/totalAtoms
        fTot=fTot/totalAtoms
        #print(MatName,sTot,pTot,dTot,fTot)
        
        Masses=[ElementData[a][4] for a in AtomicSpecies] #mass descriptors
        RedMass=(sum([b/a for a,b in zip(Masses,AtomicWeights)]))**-1
        AverageMass=np.average(Masses,weights=AtomicWeights)
        DevMass=math.sqrt(abs(np.average([c**2 for c in Masses],weights=AtomicWeights)-AverageMass**2))
        
        ENegs=[ElementData[a][5] for a in AtomicSpecies] #electronegativity descriptors
        AvgENeg=np.average(ENegs,weights=AtomicWeights)
        DiffENeg=max(ENegs)-min(ENegs)
        DevENeg=math.sqrt(abs(np.average([c**2 for c in ENegs],weights=AtomicWeights)-AvgENeg**2))
        #print(ENegs,AvgENeg,DiffENeg,DevENeg)
        
        IonEns=[ElementData[a][6] for a in AtomicSpecies]
        AtomicRads=[ElementData[a][7] for a in AtomicSpecies]
        Springs=[]
        #do not include multiplicity numerically- already taken care of by loops
        for q in range(0,len(AtomicSpecies)): #this will be higher electronegativity element
            for qq in range(0,len(AtomicSpecies)): #this will be lower electronegativity element
                if ENegs[q]>=ENegs[qq]:
                    ENHigh=ENegs[q] #electronegativity of more electronegative
                    IonEnLow=IonEns[qq] #ionization energy of less electronegative
                    R1=AtomicRads[q]
                    R2=AtomicRads[qq]
                    Springs.append((ENHigh+IonEnLow)/((R1+R2)**2))
                    #freq=math.sqrt((ENHigh+IonEnLow)/((R1+R2)**2)/RedMass)
        Omegas=[math.sqrt(a/RedMass) for a in Springs]
        OmegaAvg=np.mean(Omegas)
        OmegaDev=math.sqrt(abs(np.average([c**2 for c in Omegas])-OmegaAvg**2))
        OmegaMax=max(Omegas)
        OmegaMin=min(Omegas)
        #print(OmegaAvg,OmegaDev,OmegaMax,OmegaMin)
        
        
        IonicRadii=[ElementData[a][8] for a in AtomicSpecies]
        AvgIonicRad=np.mean(IonicRadii)
        DevIonicRad=math.sqrt(abs(np.average([c**2 for c in IonicRadii])-AvgIonicRad**2))
        
        Dipoles=[ElementData[a][9] for a in AtomicSpecies]
        AvgDipole=np.mean(Dipoles)
        DevDipole=math.sqrt(abs(np.average([c**2 for c in Dipoles])-AvgDipole**2))
        DiffDipole=max(Dipoles)-min(Dipoles)
        
        block = {'s':0,'p':1,'d':2,'f':3}
        Block=[ElementData[a][10] for a in AtomicSpecies]
        AvgBlock=np.mean(Block)
        DiffBlock=max(Block)-min(Block)
        
        Volumes=[ElementData[a][11] for a in AtomicSpecies] #neglecting prefactor
        TotalVolume=sum([a*b for a,b in zip(Volumes,AtomicWeights)]) #include scaling
        
        
        #Density of States####delete
        DOSes=[ElementData[a][12] for a in AtomicSpecies]
        TotalDOS=sum([a*b for a,b in zip(DOSes,AtomicWeights)])
        DOSperVol=TotalDOS/TotalVolume
        AvgDOS=TotalDOS/totalAtoms
        MaxDOS=max(DOSes)
        MinDOS=min(DOSes)
        
        Susceptibilities=[ElementData[a][13] for a in AtomicSpecies]
        TotalSusc=sum([a*b for a,b in zip(Susceptibilities,AtomicWeights)])
        AvgSusc=np.mean([a*b for a,b in zip(Susceptibilities,AtomicWeights)])
        AvgAtSusc=TotalSusc/totalAtoms
        AvgVolSusc=TotalSusc/TotalVolume
        MaxSusc=max(Susceptibilities)
        MinSusc=min(Susceptibilities)
        devSusc=math.sqrt(abs(np.average([c**2 for c in Susceptibilities])-AvgSusc**2))
        
        
        data=[Comp,Compounds[i],sTot,pTot,dTot,fTot,RedMass,AverageMass,DevMass,AvgENeg,DiffENeg,DevENeg,OmegaAvg,OmegaDev,OmegaMax,OmegaMin,
              AvgIonicRad,DevIonicRad,AvgDipole,DevDipole,DiffDipole,AvgBlock,DiffBlock,TotalVolume,totalAtoms,
              math.sqrt((sTot+pTot+dTot+fTot)/TotalVolume),np.sum(Masses)/TotalVolume,
              TotalDOS,DOSperVol,AvgDOS,MaxDOS,MinDOS,TotalSusc,AvgAtSusc,AvgVolSusc,MaxSusc,MinSusc,devSusc,
              gap,symm,
              volume,formE,lattA,lattB,lattC,lattAlpha,lattBeta,lattGamma,density,NAtom
              ]
        for item in data:
            f.write(str(item)+',')
        f.write('\n')
        #f.write(str(TcValues[i])+'\n')
    except:
        data=[Comp,Compounds[i],0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,0,
              0,0,
              0,0,0,0,0,0,0,0,0,0,0,
              0,'Not',
              0,0,0,0,0,0,0,0,0,0
              ]
        for item in data:
            f.write(str(item)+',')
        f.write('\n')
        continue















f.close()



print("--- %s seconds ---" % (time.time() - start_time))
