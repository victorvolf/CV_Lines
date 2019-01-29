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

t = Table(names=('Name',
				 'H_alpha Left', 'H_alpha 1', 'H_alpha 2', 'H_alpha Right', 'H_alpha v',
				 'H_beta Left',  'H_beta 1',  'H_beta 2',  'H_beta Right',  'H_beta v',
				 'HeI Left',     'HeI 1',     'HeI 2',     'HeI Right',     'HeI v',
				 'HeII Left',    'HeII 1',    'HeII 2',    'HeII Right',    'HeII v',
				 'ArIII Left',   'ArIII 1',   'ArIII 2',   'ArIII Right',   'ArIII v',
				 'class'), 
		  dtype=('U15',
                 'f4','f4','f4','f4','f4',
                 'f4','f4','f4','f4','f4',
                 'f4','f4','f4','f4','f4',
                 'f4','f4','f4','f4','f4',
		         'f4','f4','f4','f4','f4','i2'))

cv_t = fits.open('tableCV.fits', memmap = True)
cv_names = Table(cv_t[1].data)['Name']

def plot_spectra(w, f):
    nb,d = np.histogram(w,args.nbin)
    wb,d = np.histogram(w,args.nbin, weights=w)
    fb,d = np.histogram(w,args.nbin, weights=f)
    ok = nb > 0
    wb, fb = wb[ok] / nb[ok], fb[ok] / nb[ok] 
    fhi = 1.1 * np.percentile(fb, 99.9)
	
    rand = random.randint(0, 2)
    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b', color = 'purple')
    title = 'Lines of %s' %(spec[spec.rfind('/') + 1:])
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()

def flux_var(wavelength):
    _index = np.nonzero(w > wavelength)[0][0]
    range_1 = range(_index - 25, _index + 25)
    range_2 = range(_index - 5, _index + 5)
    left_range = range(_index - 35, _index - 25)
    right_range = range(_index + 25, _index + 35)
    
    if 0 in iv[range_1]:
        t_var1 = 0
    else:
        t_var1 = sum((1/iv[range_1]) * (wd[range_1] ** 2))
    t_flux1 = sum(wd[range_1] * f[range_1])
    t_flux2 = sum(wd[range_2] * f[range_2])
    t_fluxL = sum(wd[left_range] * f[left_range])
    t_fluxR = sum(wd[right_range] * f[right_range])
    
    return t_flux1, t_flux2, t_var1, t_fluxL, t_fluxR
    
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
    he2 = flux_var(4686)
    h_b = flux_var(4862)
    he1 = flux_var(5876)
    h_a = flux_var(6561)
    ar3 = flux_var(7135)
    
    t.add_row([spec[spec.find('-') + 1:spec.find('.')],
              h_a[1] - h_a[3], h_a[0], h_a[1], h_a[1] - h_a[4], h_a[2],
              h_b[1] - h_b[3], h_b[0], h_b[1], h_b[1] - h_b[4], h_b[2],
              he1[1] - he1[3], he1[0], he1[1], he1[1] - he1[4], he1[2],
              he2[1] - he2[3], he2[0], he2[1], he2[1] - he2[4], he2[2],
              ar3[1] - ar3[3], ar3[0], ar3[1], ar3[1] - ar3[4], ar3[2],
			  cls])
    '''for e in (h_a, h_b, s_2, he1, he2):
        plot_spectra(w[e[4] - 50:e[4] + 50],
                     f[e[4] - 50:e[4] + 50])'''
t.write("table.fits", "w")
print("%.3fs" % (time.time() - start))
# tyj3i5MJQseF
