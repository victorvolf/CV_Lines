#!/usr/bin/env python

usage = \
"""
Plots SDSS DR12 spectra
"""
import time
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

t = Table(names=('Name',
				 'H-alpha Left', 'H-alpha', 'H-alpha Right','H-alpha v',
				 'H-beta Left',  'H-beta',  'H-beta Right', 'H-beta v',
				 'H-gama Left',  'H-gama',  'H-gama Right', 'H-gama v',
				 'HeI Left',     'HeI',     'HeI Right',    'HeI v',
				 'HeII Left',    'HeII',    'HeII Right',   'HeII v', 
				 'class'), 
		  dtype=('U15','f4','f4','f4','f4','f4','f4','f4','f4', 'f4','f4',
		         'f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','i2'))

cv_t = fits.open('tableCV.fits', memmap = True)
cv_names = Table(cv_t[1].data)['Name']

def plot_spectra(w, f):
    nb,d = np.histogram(w,args.nbin)
    wb,d = np.histogram(w,args.nbin, weights=w)
    fb,d = np.histogram(w,args.nbin, weights=f)
    ok = nb > 0
    wb, fb = wb[ok] / nb[ok], fb[ok] / nb[ok] 
    fhi = 1.1 * np.percentile(fb, 99.9)
	
    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b')
    title = 'Lines of %s' %(spec[spec.rfind('/') + 1:])
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()

def flux_var(wavelength):
    _index = np.nonzero(w > wavelength)[0][0]
    _range = range(_index - 25, _index + 25)
    left_range = range(_index - 50, _index - 25)
    right_range = range(_index + 25, _index + 50)
    
    t_flux = sum(wd[_range] * f[_range])
    t_fluxL = sum(wd[left_range] * f[left_range])
    t_fluxR = sum(wd[right_range] * f[right_range])
    if 0 in iv[_range]:
        t_var = 0
    else:
        t_var = sum((1/iv[_range]) * (wd[_range] ** 2))
    
    return t_flux, t_var, t_fluxL, t_fluxR, _index

#w_total, f_total, ix = [0]*len(sdss), [0]*len(sdss), 0
timp = time.time()
for spec in sdss:
    hdul = fits.open(spec)
    if (all(hdul[1].data['flux']) == 0):
        continue
	
    #print(hdul[3].columns)

    table = hdul[1].data
    iv = table['ivar']
    ok = iv > 0.
    
    w  = subs.vac2air(10. ** table['loglam'])[ok]
    f  = table['flux'][ok]
    wd = table['wdisp'][ok]
    
    #w_total[ix] = w
    #f_total[ix] = f
    #ix += 1
    
    if spec[spec.find('-') + 1:spec.find('.')] in cv_names:
        cls = 1
    else:
        cls = 0
    h_a = flux_var(6561)
    h_b = flux_var(4862)
    h_g = flux_var(4341)
    he1 = flux_var(5876)
    he2 = flux_var(4686)
    
    t.add_row([spec[spec.find('-') + 1:spec.find('.')],
			  h_a[2], h_a[0], h_a[3], h_a[1],
			  h_b[2], h_b[0], h_b[3], h_b[1],
			  h_g[2], h_g[0], h_g[3], h_g[1],
			  he1[2], he1[0], he1[3], he1[1],
			  he2[2], he2[0], he2[3], he2[1],
			  cls])	   
			  
    
    '''plot_spectra(w,f)
    for e in (h_a, h_b, h_g, he1, he2):
        plot_spectra(w[e[4] - 50:e[4] + 50],
                     f[e[4] - 50:e[4] + 50])'''
    
#w_total = sorted([sum(x)/np.count_nonzero(x) for x in longest(*w_total, fillvalue = 0)])
#f_total = [sum(x)/np.count_nonzero(x) for x in longest(*f_total, fillvalue = 0)]
#plot_spectra(w_total,f_total)

print(t)
t.write("table.fits", "w")
print("%.3fs" % (time.time()-timp))
