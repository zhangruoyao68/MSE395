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
import csv
from pymatgen.ext.matproj import MPRester, MPRestError
#from sklearn import linear_model
#import matplotlib.pyplot as plt

start_time = time.time()
#3/7/18

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

#test.cv is the compound list
MaterialsFile=pd.read_csv("test_static.csv") #test1.csv is the modified data we collected
print(MaterialsFile)
#
#g=open('CompoundDescriptorsMP_exception.csv','w')
f=open('CompoundDescriptorsMP_static.csv','w')
ID=[]
target=[]

for i in range(0,len(MaterialsFile)):
  ID.append(MaterialsFile['Materials ID'][i]) #containing all the IDs of training set
  target.append(MaterialsFile['Epsilon'][i])
numMaterials=len(ID)
print("reading SourceFile completion")

#print(ID)
print("Number Of Materials Analyzed:---%s---" % (numMaterials))

API_Key='pF2RrzEPyZmIsLST'
m=MPRester(API_Key)
'''
MPID_DICT= m.query(criteria={}, properties=["task_id"])
numMaterials=len(MPID_DICT)
print("Number Of Materials Analyzed:---%s---" % (numMaterials))
#print(ID)
#print("Number Of Materials Analyzed:---%s---" % (numMaterials))
'''
Compounds=[]
for i in range(0,len(ID)):
    stabdata = m.query(criteria={'material_id':ID[i]}, properties=["pretty_formula"])
    print(stabdata[0])
    Compounds.append(stabdata[0]['pretty_formula'])
#print(Compounds)

print("API established")
ElementList=np.genfromtxt('NameListFile.txt',dtype=None)
Elements=[]
for i in range(0,len(ElementList)):
    Elements.append(ElementList[i][1].decode("utf-8"))
ElemData=[]#[[]]*len(Elements)
#print(Elements)
#print(element(8).ionic_radii)
#print("Periodic Table Completion ---%s---" %(time.time()-start_time))
block = {'s':0,'p':1,'d':2,'f':3}

for elem in Elements:
    row=[elem]
    sTot=0
    pTot=0
    dTot=0
    fTot=0
    elconfig=(element(elem).econf)
    #print(elconfig)
    sEl=re.findall(r'(?<=s)(\d*)',elconfig)
    pEl=re.findall(r'(?<=p)(\d*)',elconfig)
    dEl=re.findall(r'(?<=d)(\d*)',elconfig)
    fEl=re.findall(r'(?<=f)(\d*)',elconfig)
    #print(elconfig,sEl,pEl,dEl,fEl)

    if sEl!=[]:
        if sEl==['']:
            sTot+=1
        else:
            sTot+=int(sEl[0])
    if pEl!=[]:
        if pEl==['']:
            pTot+=1
        else:
            pTot+=int(pEl[0])
    if dEl!=[]:
        if dEl==['']:
            dTot+=1
        else:
            dTot+=int(dEl[0])
    if fEl!=[]:
        if fEl==['']:
            fTot+=1
        else:
            fTot+=int(fEl[0])
    #print(sTot,pTot,dTot,fTot)

    row.extend([sTot,pTot,dTot,fTot])
    row.append(element(elem).mass)
    row.append(element(elem).en_pauling)
    row.append(element(elem).ionenergies.get(1))

    #atomic radius for O(8), F(9), Cl(17), I(53),Br,Pm and is missing
    #manual input here
    atom_r = 0.0;
    if elem=='O':
        atom_r=65.0
    elif elem=='F':
        atom_r=57.0
    elif elem=='Cl':
        atom_r=97.0
    elif elem=='I':
        atom_r=132.0
    elif elem=='Br':
        atom_r=112.0
    elif elem=='Pm':
        atom_r=262.0
    else:
        atom_r=element(elem).atomic_radius



    row.append(atom_r)


    #!!!ionic radius data unrealiable
    row.append(element(elem).ionic_radii[0].ionic_radius)

    row.append(element(elem).dipole_polarizability)
    row.append(float(block[element(elem).block]))
    row.append(atom_r**3*10**(-7))
    ElemData.append(row)
    #print(row)
#ElementData={d[0]: d[1:] for d in ElemData}
#print(ElemData)
Database=[]
print("Element Data Collection Completion:--- %s seconds ---" % (time.time() - start_time))
    #csvfile.write(Comp+,"+sTot+","+pTot+","+dTot+","+fTor,AverageMass,DevMass,RedMass,AvgENeg,DiffENeg,DevENeg,AvgDipole')
    #csvfile.write('DevDipole,DiffDipole,TotalVolume,totalAtoms,density,OmegaAvg,OmegaDev,OmegaMax,OmegaMin')
    #csvfile.write('eHull,formE,gap,symm,lattA,lattB,lattC,lattAlpha,lattBeta,lattGamma\n')
f.write('Comp,Epsilon,sTot,pTot,dTot,fTot,NAtom,AverageMass,DevMass,RedMass,AvgENeg,DiffENeg,DevENeg,AvgDipole,DevDipole,DiffDipole,TotalVolume,totalAtoms,density,OmegaAvg,OmegaDev,OmegaMax,OmegaMin,eHull,formE,gap,symm,lattA,lattB,lattC,lattAlpha,lattBeta,lattGamma')
f.write('\n')

numberofanalyzed=0
for i in range(0,numMaterials):
    #Comp,MPID=Compounds[i].split("_")
    #Comp=Compounds[i]
    #print(Comp)
    #MatName=NameAnalyzer(Comp)
    #if 'Ac' in MatName:
    #    continue

    #res=str(MPID_DICT[i])
    #MPID_length=len(res)
    stabdata = m.query(criteria={'material_id':ID[i]}, properties=["pretty_formula", "e_above_hull","formation_energy_per_atom","band_gap","spacegroup.number",'volume',
    'structure','density','e_above_hull','nsites'])

    #print(stabdata)
    #print(stabdata[0]['spacegroup.number'])
    #stabdata = m.query(criteria={"pretty_formula": Comp}, properties=["material_id"])
    #if stabdata:
        #mini=min(stabdata, key=lambda x:x['formation_energy_per_atom'])
        #if 'material_id' in mini:
    #print(stabdata)
    symm=SymmGroupFinder(stabdata[0]['spacegroup.number'])
    #num = Compounds[i]
    gap = stabdata[0]['band_gap']
    volume=stabdata[0]['volume']
    formE=stabdata[0]['formation_energy_per_atom']
    '''
    lattA=str.format("{0:.6f}",stabdata[0]['structure'].lattice.a)
    lattB=str.format("{0:.6f}",stabdata[0]['structure'].lattice.b)
    lattC=str.format("{0:.6f}",stabdata[0]['structure'].lattice.c)
    lattAlpha=str.format("{0:.6f}",stabdata[0]['structure'].lattice.angles[0])
    lattBeta=str.format("{0:.6f}",stabdata[0]['structure'].lattice.angles[1])
    lattGamma=str.format("{0:.6f}",stabdata[0]['structure'].lattice.angles[2])
    '''
    lattA=stabdata[0]['structure'].lattice.a
    lattB=stabdata[0]['structure'].lattice.b
    lattC=stabdata[0]['structure'].lattice.c
    lattAlpha=stabdata[0]['structure'].lattice.angles[0]
    lattBeta=stabdata[0]['structure'].lattice.angles[1]
    lattGamma=stabdata[0]['structure'].lattice.angles[2]

    density=stabdata[0]['density']
    eHull=stabdata[0]['e_above_hull']
    NAtom=stabdata[0]['nsites']
    Comp=stabdata[0]['pretty_formula']
    MatName=NameAnalyzer(Comp)
    #if ('Ac' or 'Pu') in MatName:
    #    continue
    #print(Comp,num,gap)

    try:
        #print("try")
        AtomicRatio=[float(k) for k in MatName[1::2]]
        #print(AtomicRatio)
        AtomicSpecies=[k for k in MatName[0::2]]
        totalAtoms=sum([float(k) for k in MatName[1::2]]) #total number of atoms in formula unit
        #print(MatName,AtomicSpecies,AtomicWeights,totalAtoms)

        #initization of all data to be pulled from ElemData
        sTot=0
        pTot=0
        dTot=0
        fTot=0
        Masses=[]
        ENeg=[]
        IonEns=[]
        AtomicRads=[]
        Dipoles=[]
        Volumes=[]

        #random initization for highest and lowest E
        ENHigh=0
        IonEnLow=100
        AtomicRads_ENHigh=0
        AtomicRads_IonEnLow=0

        
        
        for j in range(0,len(AtomicSpecies)): #this gets reduced electronic configuration
            #print(AtomicSpecies[j])

            for x in range(0,len(ElemData)):

                if ElemData[x][0]==AtomicSpecies[j]:
                    
                    sTot+=ElemData[x][1]*AtomicRatio[j]
                    pTot+=ElemData[x][2]*AtomicRatio[j]
                    dTot+=ElemData[x][3]*AtomicRatio[j]
                    fTot+=ElemData[x][4]*AtomicRatio[j]
                    for y in range(0,int(AtomicRatio[j])):
                        Masses.append(ElemData[x][5])
                        ENeg.append(ElemData[x][6])
                        IonEns.append(ElemData[x][7])
                        AtomicRads.append(ElemData[x][8])
                        Dipoles.append(ElemData[x][10])
                        Volumes.append(ElemData[x][12])
                    #if ElemData[x][6]>ENHigh:
                        #ENHigh=ElemData[x][6]
                        #AtomicRads_ENHigh=ElemData[x][8]
                    #if ElemData[x][6]<IonEnLow:
                        #IonEnLow=ElemData[x][6]
                        #AtomicRads_IonEnLow=ElementData[x][8]
                        '''
        if (find_all<len(AtomicSpecies)):
            g.write(Comp)
            g.write('\n')
            '''


        row1=[]

        sTot=sTot/totalAtoms

        pTot=pTot/totalAtoms
        #row1.append(pTot)
        dTot=dTot/totalAtoms
        #row1.append(dTot)
        fTot=fTot/totalAtoms
        #row1.append(fTot)

        AverageMass=np.average(Masses)
        DevMass=np.std(Masses)
        RedMass=1/(np.sum(1./np.array(Masses)))

        AvgENeg=np.average(ENeg)
        DiffENeg=max(ENeg)-min(ENeg)
        DevENeg=np.std(ENeg)

        AvgDipole=np.average(Dipoles)
        DevDipole=np.std(Dipoles)
        DiffDipole=max(Dipoles)-min(Dipoles)
        TotalVolume = np.sum(Volumes)

        #print("Halfway")
        '''
        Masses=[ElemData[a][5] for a in AtomicSpecies] #mass descriptors
        #RedMass=(sum([b/a for a,b in zip(Masses,AtomicRatio)]))**-1
        #AverageMass=np.average(Masses,weights=AtomicRatio)
        #DevMass=math.sqrt(abs(np.average([c**2 for c in Masses],weights=AtomicRatio)-AverageMass**2))

        #ENegs=[Elem[a][5] for a in AtomicSpecies] #electronegativity descriptors
        #AvgENeg=np.average(ENegs,weights=AtomicRatio)
        #DiffENeg=max(ENegs)-min(ENegs)
        #DevENeg=math.sqrt(abs(np.average([c**2 for c in ENegs],weights=AtomicRatio)-AvgENeg**2))
        #print(ENegs,AvgENeg,DiffENeg,DevENeg)

        #IonEns=[ElementData[a][6] for a in AtomicSpecies]
        #AtomicRads=[ElementData[a][7] for a in AtomicSpecies]
        '''
        #spring constant calculation
        Springs=[]
        #do not include multiplicity numerically- already taken care of by loops
        #print(ENeg)
        for q in range(0,len(ENeg)): #this will be higher electronegativity element
            for qq in range(0,len(ENeg)): #this will be lower electronegativity element
                if ENeg[q]>=ENeg[qq]:
                    ENHigh=ENeg[q] #electronegativity of more electronegative
                    IonEnLow=ENeg[qq] #ionization energy of less electronegative
                    AtomicRads_ENHigh=AtomicRads[q]
                    AtomicRads_IonEnLow=AtomicRads[qq]
                    Springs.append((ENHigh+IonEnLow)/((AtomicRads_ENHigh+AtomicRads_IonEnLow)**2))
                    #freq=math.sqrt((ENHigh+IonEnLow)/((R1+R2)**2)/RedMass)


        #Omega calculation
        Omegas=[]
        for h in range(0,len(Springs)):
            Omegas.append(math.sqrt(Springs[h]/RedMass))

        OmegaAvg=np.mean(Omegas)
        OmegaDev=np.std(Omegas)
        OmegaMax=max(Omegas)
        OmegaMin=min(Omegas)


        #output list
        #print(Comp)
        #print(sTot,pTot,dTot,fTot)
        #print(AverageMass,DevMass,RedMass)
        #print(AvgENeg,DiffENeg,DevENeg)
        #print(AvgDipole,DevDipole,DiffDipole)
        #print(TotalVolume,totalAtoms,density)
        #print(OmegaAvg,OmegaDev,OmegaMax,OmegaMin)
        #print(eHull,formE,gap)
        #print(symm,lattA,lattB,lattC,lattAlpha,lattBeta,lattGamma)



        row=[Comp,target[i],sTot,pTot,dTot,fTot,NAtom,AverageMass,DevMass,RedMass,AvgENeg,DiffENeg,DevENeg,AvgDipole,DevDipole,DiffDipole,TotalVolume,totalAtoms,density,OmegaAvg,OmegaDev,OmegaMax,OmegaMin,eHull,formE,gap,symm,lattA,lattB,lattC,lattAlpha,lattBeta,lattGamma]
        #print(row1)
        for item in row:
            f.write(str(item)+',')
        f.write('\n')
        Database.append(row)
        print("Compound Analyzed:"+Comp)
        numberofanalyzed+=1





        #ion radius data unrealiable
        '''
        IonicRadii=[ElementData[a][8] for a in AtomicSpecies]
        AvgIonicRad=np.mean(IonicRadii)
        DevIonicRad=math.sqrt(abs(np.average([c**2 for c in IonicRadii])-AvgIonicRad**2))


        Dipoles=[ElementData[a][9] for a in AtomicSpecies]
        AvgDipole=np.mean(Dipoles)
        DevDipole=math.sqrt(abs(np.average([c**2 for c in Dipoles])-AvgDipole**2))
        DiffDipole=max(Dipoles)-min(Dipoles)
        print(AvgDipole)


        block = {'s':0,'p':1,'d':2,'f':3}
        Block=[ElementData[a][10] for a in AtomicSpecies]
        AvgBlock=np.mean(Block)
        DiffBlock=max(Block)-min(Block)
        '''

        #Volumes=[ElementData[a][11] for a in AtomicSpecies] #neglecting prefactor
        #TotalVolume=sum([a*b for a,b in zip(Volumes,AtomicRatio)]) #include scaling


        #Density of States####delete
        #have not realized yet
        '''
        DOSes=[ElementData[a][12] for a in AtomicSpecies]
        TotalDOS=sum([a*b for a,b in zip(DOSes,AtomicRatio)])
        DOSperVol=TotalDOS/TotalVolume
        AvgDOS=TotalDOS/totalAtoms
        MaxDOS=max(DOSes)
        MinDOS=min(DOSes)

        Susceptibilities=[ElementData[a][13] for a in AtomicSpecies]
        TotalSusc=sum([a*b for a,b in zip(Susceptibilities,AtomicRatio)])
        AvgSusc=np.mean([a*b for a,b in zip(Susceptibilities,AtomicRatio)])
        AvgAtSusc=TotalSusc/totalAtoms
        AvgVolSusc=TotalSusc/TotalVolume
        MaxSusc=max(Susceptibilities)
        MinSusc=min(Susceptibilities)
        devSusc=math.sqrt(abs(np.average([c**2 for c in Susceptibilities])-AvgSusc**2))

        
        data=[Comp,Compounds[i],sTot,pTot,dTot,fTot,RedMass,AverageMass,DevMass,AvgENeg,DiffENeg,DevENeg,OmegaAvg,OmegaDev,OmegaMax,OmegaMin,
              AvgIonicRad,DevIonicRad,AvgDipole,DevDipole,DiffDipole,AvgBlock,DiffBlock,TotalVolume,totalAtoms,
              math.sqrt((sTot+pTot+dTot+fTot)/TotalVolume),np.sum(Masses)/TotalVolume,
              TotalDOS,DOSperVol,AvgDOS,MaxDOS,MinDOS,TotalSusc,AvgAtSusc,AvgVolSusc,MaxSusc,MinSusc,devSusc,
              gap,symm,volume,formE,lattA,lattB,lattC,lattAlpha,lattBeta,lattGamma,density,NAtom]
        print('data')
        
        for item in data:
            f.write(str(item)+',')
        f.write('\n')
        #f.write(str(TcValues[i])+'\n')
        '''
        
    except:
        print("exception")
        data=[Comp,Comp,0,0,0,0,0,0,0,0,0,0,0,0,0,
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

#g.close()

#ar=np.asarray(Database)
#print(Database)
#ar=str(np.asarray(Database))
#np.savetxt('descriptors.csv',ar,delimiter=',')
f.close()

print("Number Of Materials Analyzed:---%s---" % (numberofanalyzed))
print("Compound Data Sourcing:--- %s seconds ---" % (time.time() - start_time))

#fitting model
#regr=linear_model.Ridge(alpha=.1)
#regr.fit(Database,target)
