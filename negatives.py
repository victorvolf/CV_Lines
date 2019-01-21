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

negatives = pd.read_csv('negatives.csv', sep = " ", names = ['name', 'class'], header = None)
def plot_spectra(w, f, target):
    nb,d = np.histogram(w,args.nbin)
    wb,d = np.histogram(w,args.nbin, weights=w)
    fb,d = np.histogram(w,args.nbin, weights=f)
    ok = nb > 0
    wb, fb = wb[ok] / nb[ok], fb[ok] / nb[ok] 
    fhi = 1.1 * np.percentile(fb, 99.9)
	
    rand = random.randint(0, 2)
    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b', color = 'blue')
    if target.values == 0:
        title = 'Lines of non-CV %s' %(spec[spec.rfind('/') + 1:])
    else:
        title = 'Lines of CV %s' %(spec[spec.rfind('/') + 1:])
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()

for spec in sdss:
    if not spec[spec.find('-') + 1:spec.find('.')] in list(negatives['name']):
        continue
    else:
        hdul = fits.open(spec)
        if (all(hdul[1].data['flux']) == 0):
            continue
        table = hdul[1].data
        iv = table['ivar']
        ok = iv > 0.
        
        w  = subs.vac2air(10. ** table['loglam'])[ok]
        f  = table['flux'][ok]
        plot_spectra(w,f,negatives[negatives['name'] == spec[spec.find('-') + 1:spec.find('.')]]['class'])

print("%.3fs" % (time.time() - start))
