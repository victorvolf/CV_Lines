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

from sklearn import ensemble as ens
from sklearn import model_selection as cv
from sklearn.metrics import roc_curve, auc
from sklearn.utils import shuffle

def perf_measure(y_actual, y_hat):
    TP, FP, TN, FN = 0, 0, 0, 0

    for i in range(len(y_hat)): 
        if y_actual[i] == y_hat[i] == 1:
           TP += 1
        if y_hat[i] == 1 and y_actual[i] != y_hat[i]:
           FP += 1
        if y_actual[i] == y_hat[i] == 0:
           TN += 1
        if y_hat[i] == 0 and y_actual[i] != y_hat[i]:
           FN += 1
    return TP, FP, TN, FN

t = fits.open('table.fits', memmap = True)
table = Table(t[1].data)

df = table.to_pandas()
df = shuffle(df)

classT = df['class'].to_frame()
df = df.drop(['Name'], axis=1)
df = df.drop(['class'], axis=1)

kf = cv.KFold(10)
avg_acc = 0
imp = [0] * 25
performance=[0, 0, 0, 0]

for train_index, test_index in kf.split(df):
    xTrain, xTest = df.iloc[train_index, :], df.iloc[test_index, :]
    yTrain, yTest = classT.iloc[train_index, :], classT.iloc[test_index, :]
    
    rfc = ens.RandomForestClassifier(n_estimators = 100,
                                     n_jobs = -1, 
                                     criterion = 'entropy',
                                     max_features = 'sqrt')
    
    rfc.fit(xTrain.values, yTrain.values.ravel())
    predictions = rfc.predict(xTest)
    fpr,tpr,threshold = roc_curve(yTest['class'], predictions, pos_label = 1)
    
    p= perf_measure(list(yTest['class']), predictions)
    
    for i in range(0, 4):
        performance[i] += p[i]
    
    imp += rfc.feature_importances_
    avg_acc += auc(fpr, tpr)
avg_acc /= 10

print(' TP   FP  TN    FN')
print(performance)
print(avg_acc)
print(imp)
