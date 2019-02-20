#!/usr/bin/env python

usage = \
"""
Plots SDSS DR12 spectra
"""
import time
import os, sys, argparse
import numpy as np
import pandas as pd
import random
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import astropy.io.fits as fits
from astropy.table import Table, Column
from trm import subs
from itertools import zip_longest as longest

parser = argparse.ArgumentParser(description = usage)

# positional
parser.add_argument('sdss', 
                     help = 'name of SDSS FITS file(s)', 
                     nargs = "*")

# optional
parser.add_argument('-n', 
                    dest = 'nbin', 
                    type = int,
                    default = 1000, 
                    help = 'number of pixels to bin the data into')

# OK, done with arguments.
args = parser.parse_args()
sdss = args.sdss
start = time.time()

n = 100

t = Table([[1]])
t['Name'] = "abcdefghijklmno"
for i in range(0, n):
    t['%sf' %(i)] = 0.
    t['%sv' %(i)] = 0.
t['class'] = 0
t.remove_column('col0')
t.remove_row(0)

cv_t = fits.open('tableCV.fits', memmap = True)
cv_names = Table(cv_t[1].data)['Name']

def flux_var(n):
    input_space = [spec[spec.find('-') + 1:spec.find('.')]]
    
    min_index = np.nonzero(w > min(4500, w[-1]))[0][0]
    max_index = np.nonzero(w > max(7000, w[-2]))[0][0]
    
    chunks_var = np.array_split(iv[min_index:max_index], n)
    chunks_wd = np.array_split(wd[min_index:max_index], n)
    chunks_f = np.array_split(f[min_index:max_index], n)
    #print(len(chunks_var),len(chunks_wd),len(chunks_f))
    
    for i in range(0, n):
        if 0 in chunks_var[i]:
            t_var = 0
        else:
            t_var = sum((1/chunks_var[i]) * (chunks_wd[i] ** 2))
        t_flux = sum(chunks_wd[i] * chunks_f[i])    
        input_space.extend([t_flux, t_var])
    input_space.append(cls)
    t.add_row(input_space)
   
for spec in sdss:
    hdul = fits.open(spec)
    if all(hdul[1].data['flux']) == 0:
        continue
        
    table = hdul[1].data
    iv = table['ivar']
    ok = iv > 0.
    
    w  = subs.vac2air(10. ** table['loglam'])[ok]
    f  = table['flux'][ok]
    wd = table['wdisp'][ok]
    if spec[spec.find('-') + 1:spec.find('.')] in cv_names:
        cls = 1
    else:
        cls = 0
    
    flux_var(n)
print(t)
t.write("table2.fits", "w")
print("%.3fs" % (time.time() - start))
# tyj3i5MJQseF
