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
from sklearn import linear_model
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.neural_network import MLPRegressor as ML
from sklearn import metrics

def breakingdown(Matrix,iteration,skipSize):
    new_Matrix=[]
    for item in range(0,len(Matrix)):
        if item<iteration*skipSize or item>(iteration+1)*skipSize:
            new_Matrix=np.append(Matrix(item))
    return new_Matrix


start_time=time.time()

MaterialFile=pd.read_csv("CompoundDescriptorsMP_static.csv")
f=open('NeuralNetworkRegressorResults_static.csv','w')
#MaterialFile2= pd.read_csv("CoulombResults.csv")
#print("Input File reading Complete")
target=[]
descriptors=[]
Name=[]
target_test=[]
descriptors_test=[]
TrainingName=[]
TestingName=[]
for i in range(0,len(MaterialFile)):
    #if(MaterialFile['LattAlpha'][i]==MaterialFile['LattBeta'][i]==MaterialFile['LattGamma'][i]=90)
    row1=[]
    #Name.append(MaterialFile2['Material ID'][i])
    '''
    should include CM later
    #row1.append(MaterialFile2['cm0'][i]/1000)
    #row1.append(MaterialFile2['cm1'][i]/1000)
    #row1.append(MaterialFile2['cm2'][i]/1000)
    '''
    row1.append(MaterialFile['AvgIonicRad'][i])
    row1.append(MaterialFile['DevIonicRad'][i])
    row1.append(MaterialFile['VolumeAverageMass'][i])
    row1.append(MaterialFile['AvgElectron'][i])
    row1.append(MaterialFile['NAtom'][i])
    row1.append(MaterialFile['sTot'][i])
    row1.append(MaterialFile['pTot'][i])
    row1.append(MaterialFile['dTot'][i])
    row1.append(MaterialFile['fTot'][i])
    row1.append(MaterialFile['AverageMass'][i])
    row1.append(MaterialFile['DevMass'][i])
    row1.append(MaterialFile['RedMass'][i])
    row1.append(MaterialFile['AvgENeg'][i])
    row1.append(MaterialFile['DevENeg'][i])
    row1.append(MaterialFile['DiffENeg'][i])
    row1.append(MaterialFile['AvgDipole'][i])
    row1.append(MaterialFile['DevDipole'][i])
    row1.append(MaterialFile['DiffDipole'][i])
    row1.append(MaterialFile['TotalVolume'][i])
    row1.append(MaterialFile['totalAtoms'][i])
    row1.append(MaterialFile['density'][i])
    row1.append(MaterialFile['OmegaAvg'][i])
    row1.append(MaterialFile['OmegaDev'][i])
    row1.append(MaterialFile['gap'][i])
    row1.append(MaterialFile['formE'][i])
    row1.append(MaterialFile['lattA'][i])
    row1.append(MaterialFile['lattB'][i])
    row1.append(MaterialFile['lattC'][i])
    row1.append(MaterialFile['lattAlpha'][i])
    row1.append(MaterialFile['lattBeta'][i])
    row1.append(MaterialFile['lattGamma'][i])
    if MaterialFile['lattAlpha'][i]>0: #==MaterialFile['lattBeta'][i]==MaterialFile['lattGamma'][i]:
        if i%10!=0:
            descriptors.append(row1)
            target.append(MaterialFile['Epsilon'][i])
            TrainingName.append(MaterialFile['Comp'][i])
        else:
            descriptors_test.append(row1)
            target_test.append(MaterialFile['Epsilon'][i])
            TestingName.append(MaterialFile['Comp'][i])
print("Input File reading Complete")
print(len(descriptors))

#chossing the best alpha

#regr=linear_model.Lasso()
#regr=linear_model.LinearRegression()
#alphas=[0.001,0.01,0.1,1,10,12,100,1000]
#scores=[regr.set_params(alpha=alpha).fit(descriptors,target).score(descriptors_test,target_test) for alpha in alphas]
#print(scores)
#best_alpha = alphas[scores.index(max(scores))]
#regr.alpha = best_alpha
#fit=regr.fit(descriptors, target)

#regr=RandomForestRegressor(n_estimators=,criterion='mse',warm_start=True)

#cross-validation
#for iteration in range(0,10):
    #if iteration==0:
        #first_regr=regr.fit(breakingdown(descriptors,iteration,10),breakingdown(target,iteration,10))
    #else:
        #re_regr=RandomForestRegressor(warm_start=True)
        #fitting=reregr.fit(breakingdown(descriptors,iteration,10),breakingdown(target,iteration,10))
iteration_cycle=10
iteration=0
n_estimators=15
#kf=KFold(n_split=iteration_cycle)
regr=RandomForestRegressor(n_estimators=n_estimators,criterion='mse',warm_start=False)
#regr=ML()
'''

while iteration < iteration_cycle:
    n_estimators+=1
    #regr=RandomForestRegressor(n_estimators=n_estimators,criterion='mse',warm_start=True)
    x_training,x_cross,y_training,y_cross=train_test_split(descriptors,target,test_size=0.2,random_state=42)
    fitting=regr.fit(x_training,y_training)
    iteration+=1
'''
#scores=cross_val_score(regr,descriptors,target, cv=10)
fitting=regr.fit(descriptors,target)

#print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))

fitting_score=fitting.score(descriptors_test,target_test)
print(fitting_score)
#print(regr.coef_)
std_score=0
predicted_descriptor=[]
predicted_test=[]
predicted_descriptor=fitting.predict(descriptors)
predicted_test=fitting.predict(descriptors_test)
#predicted_value.append(fitting.predict(descriptors_test))
#print(predicted_descriptor)
#print(predicted_test)

for i in range(0,len(predicted_descriptor)):

    #for item in range(0,len(descriptors[i])):
    #print(predicted_value)
    #predicted_value=fitting.predict(map(list,zip(*descriptors[i])))
    f.write(str(predicted_descriptor[i])+","+str(target[i])+"\n")
    std_score+=(predicted_descriptor[i]-target[i])**2

for i in range(0,len(descriptors_test)):
    #predicted_value=fitting.predict(map(list,zip(*descriptors_test[i])))
    #for item in range(0,len(descriptors_test[i])):
        #predicted_value+=descriptors_test[i][item]*regr.coef_[item]
    print(predicted_test[i])
    print(target_test[i])
    print('\n')
    f.write(str(predicted_test[i])+","+str(target_test[i])+"\n")
    std_score+=(predicted_test[i]-target_test[i])**2

#f.write(str(fitting.predict(descriptors)))
#f.write(str(fitting.predict(descriptors_test)))
print(std_score)

#plotting
plt.plot(target,predicted_descriptor,'g^',target_test,predicted_test,'bs')
plt.xlabel('experimental value')
plt.ylabel('predicted value')
#regreession_line=np.polyfit(target,predicted_descriptor,1)
#pp=np.poly1d(regreession_line)
#plt.plot(target,pp(target),'r--')
plt.axis([0,400,0,400])
plt.show()

f.close()
