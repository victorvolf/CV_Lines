#!/usr/bin/env python
import os, sys, argparse
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import astropy.io.fits as fits
from astropy.table import Table, Column
from trm import subs
from itertools import zip_longest as longest
from collections import defaultdict

from PRF import prf
from sklearn import ensemble as ens
from sklearn import model_selection as cv
from sklearn.metrics import roc_curve, auc
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler

def perf_measure(y_actual, y_hat, nTest):
    TP, FP, TN, FN = 0, 0, 0, 0
    names_fp = []
    names_fn = []

    for i in range(len(y_hat)): 
        if y_actual[i] == y_hat[i] == 1:
            TP += 1
        if y_hat[i] == 1 and y_actual[i] != y_hat[i]:
            names_fp.append(nTest.iloc[i,:].values[0])
            FP += 1
        if y_actual[i] == y_hat[i] == 0:
            TN += 1
        if y_hat[i] == 0 and y_actual[i] != y_hat[i]:
            names_fn.append(nTest.iloc[i,:].values[0])
            FN += 1
    return TP, FP, TN, FN, names_fp, names_fn

t = fits.open('table.fits', memmap = True)
table = Table(t[1].data)

df = table.to_pandas()
df = shuffle(df)

classT = df['class'].to_frame()
names = df['Name'].to_frame()
df = df.drop(['Name', 'class', 'H_beta v'], axis=1)

kf = cv.KFold(10)
avg_acc = 0
importances = defaultdict(float)
performance=[0, 0, 0, 0]
names_fp = []
names_fn = []

for train_index, test_index in kf.split(df):
    xTrain, xTest = df.iloc[train_index, :], df.iloc[test_index, :]
    yTrain, yTest = classT.iloc[train_index, :], classT.iloc[test_index, :]
    nTrain, nTest = names.iloc[train_index, :], names.iloc[test_index, :]
    
    #rfc = prf(n_estimators=10, bootstrap=True, keep_proba=0.05)
    rfc = ens.RandomForestClassifier(n_estimators = 500,
                                     n_jobs = -1, 
                                     criterion = 'entropy',
                                     max_features = None,
                                     class_weight = 'balanced_subsample',
                                     oob_score = True,
                                     random_state = 5)
    
    rfc.fit(X = xTrain.values, y = yTrain.values.ravel())
    predictions = rfc.predict(xTest)
    fpr,tpr,threshold = roc_curve(yTest['class'], predictions, pos_label = 1)
    
    p = perf_measure(list(yTest['class']), predictions, nTest)
    for i in range(0, 4):
        performance[i] += p[i]
    for name in p[4]:		
	    names_fp.append(name)
    for name in p[5]:		
	    names_fn.append(name)
    
    for feature, importance in zip(df.columns, rfc.feature_importances_):
        importances[feature] += importance
    print(auc(fpr, tpr))
    avg_acc += auc(fpr, tpr)
print(' TP   FP  TN    FN\n', performance)
print(avg_acc / 10)
for w in sorted(importances, key=importances.get, reverse=True):
    print("%s:" % w, importances[w])
    
with open('negatives.csv', "w") as outfile:
    for entries in set(names_fp):
        outfile.write(entries + " 0")
        outfile.write("\n")
    for entries in set(names_fn):
        outfile.write(entries + " 1")
        outfile.write("\n")
