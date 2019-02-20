from __future__ import division

import os, sys, gc, re
import warnings
import pickle
import time

import numpy as np
import pandas as pd
import tensorflow as tf
import keras.backend as K
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.table import Table, Column

from collections import Counter
from tensorflow.python.client import timeline
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.metrics import r2_score,roc_curve,auc, log_loss

from keras import regularizers
from keras.optimizers import Adam
from keras.models import Sequential, Model
from keras.callbacks import EarlyStopping, TensorBoard
from keras.constraints import maxnorm, nonneg, unitnorm
from keras.callbacks import ReduceLROnPlateau,ModelCheckpoint
from keras.layers import GRU, Activation, Dropout, Input, BatchNormalization, AtrousConvolution1D
from keras.layers import Dense, Flatten, Conv2D, MaxPooling2D, Conv1D, MaxPooling1D, Lambda

def to_categorical(y, num_classes=None, dtype='float32'):
    Y = []
    for i in range(0, len(y)):
        zeros = [0] * (max(y) + 1)
        zeros[y[i]] = 1
        Y.append(zeros)
    return np.array(Y)

def build_model():
    input = Input(shape=(X_train.shape[1], 6), dtype='float32', name='input0')
    output = Conv1D(256,
                    80,
                    strides = 100,
                    padding = "same",
                    kernel_initializer = 'glorot_uniform',
                    kernel_regularizer = regularizers.l2(l = 0.0001),
                    bias_regularizer = regularizers.l2(l = 0.01),
                    use_bias = True)(input)
    output = BatchNormalization()(output)
    output = Activation('tanh')(output)
    output = MaxPooling1D(pool_size = 4,
                          strides = None,
                          padding = 'same')(output)
    output = Lambda(lambda x: K.mean(x, axis=1))(output)
    output = Dense(len(classes),
                   activation = 'softmax',
                   kernel_initializer = 'normal',
                   kernel_regularizer = regularizers.l2(l=0.01),
                   bias_regularizer = regularizers.l2(l=0.01),
                   use_bias = True)(output)
    model = Model(inputs = input, outputs = output)
    return model

t = fits.open('table1.fits', memmap = True)
train = Table(t[1].data).to_pandas()
train = shuffle(train)
classes = train['class'].to_frame()
train = train.drop(['class', 'Name'], axis=1)

metadata = pd.read_csv('metadata.csv')
print('Read training set and metadata...')

num_columns = len(train.columns)
X_train = np.array([train.values])
train = train.reset_index()
class_map = dict()
for i, val in enumerate(list(classes['class'].values)):
    class_map[i] = val
targets = classes
target_map = np.zeros((targets.shape[0],))
target_map = np.array([class_map[i] for i, val in enumerate(list(classes['class'].values))])
Y = to_categorical(target_map)
y_count = Counter(target_map)
wtable = np.zeros((len(classes),))
for i in range(len(classes)):
    wtable[i] = y_count[i] / target_map.shape[0]

y_map = target_map
y_categorical = Y
folds = KFold(n_splits=10)
start = time.time()
clfs = []
oof_preds = np.zeros((len(X_train), len(classes)))

for trn_, val_ in folds.split(X_train[0]):
    x_train, y_train = X_train[trn_], Y[trn_]
    x_valid, y_valid = X_train[val_], Y[val_]

    model = build_model()
    optimizer = Adam(lr=0.001, beta_1=0.9, beta_2=0.999)
    stopping = EarlyStopping(monitor='val_loss', patience=60, verbose=0, mode='auto')

    model.compile(loss='categorical_crossentropy', optimizer=optimizer, metrics=['accuracy'])
    history = model.fit(x_train, y_train,
                        validation_data = [x_valid, y_valid],
                        epochs = 100,
                        batch_size = 512,
                        shuffle = False,
                        verbose = 1,
                        callbacks = [stopping])

    print('Loading Best Noodle')
    # Get predicted probabilities for each class
    oof_preds[val_, :] = model.predict(x_valid,batch_size=batch_size)
    clfs.append(model)

elapsed_time = time.time() - start
print("elapsed_time:", elapsed_time)
