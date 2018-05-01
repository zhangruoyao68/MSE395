import pandas as pd
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.core import periodic_table
import numpy as np
import pymatgen as mg
from scipy.linalg import eig
import csv

API_Key='pF2RrzEPyZmIsLST'
MP_all=[]

f=open("All_ID.csv",'w')
with MPRester(API_Key) as mp:
    data = mp.query(criteria={}, properties=["task_id"])
    for i in range(len(data)):
        f.write(data[i]["task_id"])
        f.write('\n')

f.close()
print(len(MP_all))
