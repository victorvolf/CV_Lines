#!/usr/bin/env python

usage = \
"""
Plots SDSS DR12 spectra
"""

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
				 'H-beta Left', 'H-beta', 'H-beta Right','H-beta v',
				 'HeI Left', 'HeI', 'HeI Right','HeI v',
				 'HeII Left', 'HeII', 'HeII Right''HeII','HeII v', 
				 'class'), 
		  dtype=('U15','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','i4'))

cv_t = fits.open('tableCV.fits', memmap = True)
cv_names = Table(cv_t[1].data)['Name']

def plot_spectra(w, f, element):
    nb,d = np.histogram(w,args.nbin)
    wb,d = np.histogram(w,args.nbin, weights=w)
    fb,d = np.histogram(w,args.nbin, weights=f)
    ok = nb > 0
    wb, fb = wb[ok] / nb[ok], fb[ok] / nb[ok] 
    fhi = 1.1 * np.percentile(fb, 99.9)
	
    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b')
    title = '%s lines of ' %element + spec[spec.rfind('/') + 1:] + ', ' + pos
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()

def flux_var(wavelength):
    _index = np.nonzero(w > wavelength)[0][0]
    _range = range(_index - 23, _index + 23)
    ok = iv[_range] > 0
    
    t_flux = sum(w_d[_range] * f[_range])
    t_var = sum((1/iv[_range][ok]) * (w_d[_range][ok] ** 2))
    
    return t_flux, t_var, _index

w_total, f_total, ix = [0]*len(sdss), [0]*len(sdss), 0

for spec in sdss:
    hdul = fits.open(spec)
    head  = hdul[0].header
    ra  = head['PLUG_RA']/15.
    dec = head['PLUG_DEC']
    pos = subs.d2hms(ra, sep = ' ', dp = 2) + ' ' + \
		  subs.d2hms(dec, sep = ' ', dp = 1, sign = True)
    #print('Position = ',pos)
    if sum(hdul[1].data['flux']) == 0:
		continue
		
   	# Store fluxes and wavelengths in arrays 
    table = hdul[1].data
    w  = subs.vac2air(10. ** table['loglam'])
    f  = table['flux']
    iv = table['ivar']
    w_d = table['wdisp']
    ok = iv > 0.
    w  = w[ok]
    f  = f[ok]
    w_d = w_d[ok]
    #w_total[ix] = w
    #f_total[ix] = f
    #ix += 1
    
    
    if spec[spec.find('-') + 1:spec.find('.')] in cv_names:
        cls = 1
    else:
        cls = 0

    ha_left, h_a, ha_right = flux_var(6511), flux_var(6561), flux_var(6611)
    hb_left, h_b, hb_right = flux_var(4812), flux_var(4862), flux_var(4912)
    he1_left, he1, he1_right = flux_var(5826), flux_var(5876), flux_var(5926)
    he2_left, he2, he2_right = flux_var(4636), flux_var(4686), flux_var(4736)
    
    t.add_row([spec[spec.find('-') + 1:spec.find('.')],
			  ha_left[0],h_a[0], ha_right[0], h_a[1],
			  hb_left[0], h_b[0], hb_right[0], h_b[1],
			  he1_left[0] ,he1[0], he1_right[0], he1[1],
			  he2_left[0] ,he2[0], he2_right[0], he2[1],
			  cls])
			   
    '''# Entire spectra
    plot_spectra(w,f, '')
    
    # H-alpha
    plot_spectra(w[h_a[2] - 50:h_a[2] + 50],
                 f[h_a[2] - 50:h_a[2] + 50],
                 'H-alpha')
    
    # H-beta 
    plot_spectra(w[h_b[2] - 50:h_b[2] + 50],
                 f[h_b[2] - 50:h_b[2] + 50],
                 'H-beta')
                 
    # HeI
    plot_spectra(w[he1[2] - 50:he1[2] + 50],
                 f[he1[2] - 50:he1[2] + 50],
                 'HeI')
                 
    # HeII
    plot_spectra(w[he2[2] - 50:he2[2] + 50],
                 f[he2[2] - 50:he2[2] + 50],
                 'HeII')'''
#w_total = sorted([sum(x)/np.count_nonzero(x) for x in longest(*w_total, fillvalue = 0)])
#f_total = [sum(x)/np.count_nonzero(x) for x in longest(*f_total, fillvalue = 0)]

#print(t)
t.write("table.fits", "w")
