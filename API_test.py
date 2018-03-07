import pandas as pd
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.core import periodic_table
import numpy as np
import pymatgen as mg
from scipy.linalg import eig

testdata=pd.read_csv('test.csv')
print(testdata['Materials ID'][0])
Compounds=[]
for i in range(0,len(testdata)):
  Compounds.append(testdata['Materials ID'][i])

print(Compounds)
API_Key='pF2RrzEPyZmIsLST'
m=MPRester(API_Key)

for i in range(0,len(testdata)):
    #Comp,MPID=Compounds[i].split("_")
    #Comp=Compounds[i]
    #print(Comp)
    #MatName=NameAnalyzer(Comp)
    #if 'Ac' in MatName:
    #   continue
    stabdata = m.query(criteria={'material_id':Compounds[i]}, properties=["pretty_formula", "e_above_hull",
                                                                      "formation_energy_per_atom","band_gap",
                                                                      "spacegroup.number",'volume',
                                                                      'structure','density','e_above_hull','nsites'])
print(stabdata)
Comp=stabdata[0]['pretty_formula']
print(Comp)