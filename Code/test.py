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
'''
API_Key='pF2RrzEPyZmIsLST'
m=MPRester(API_Key)
data = m.query(criteria={}, properties=["task_id"])
res=str(data[1])
len=len(res)
print(res[13:len-2])

stabdata = m.query(criteria={'material_id':res[13:len-2]}, properties=["pretty_formula", "e_above_hull",
                                                                      "formation_energy_per_atom","band_gap",
                                                                      "spacegroup.number",'volume',
                                                                      'structure','density','e_above_hull','nsites'])

print(stabdata)
'''
