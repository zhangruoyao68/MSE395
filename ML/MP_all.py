import pandas as pd
from pymatgen.ext.matproj import MPRester, MPRestError
from pymatgen.core import periodic_table
import numpy as np
import pymatgen as mg
from scipy.linalg import eig
import csv

API_Key='pF2RrzEPyZmIsLST'

with MPRester(API_Key) as mp:
	data = mp.query(criteria={}, properties=["task_id"])
	print(len(data))