#!/usr/bin/env python

usage = \
"""
Plots SDSS DR12 spectra
"""

import os, sys, argparse
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import astropy.io.fits as fits
from astropy.table import Table, Column
from trm import subs

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

avg_flux = 0
t = Table(names=('H-alpha', 'H-alpha σ^2', 'H-beta', 'H-beta σ^2', 'HeII', 'HeII σ^2'), 
		  dtype=('f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

def plot_spectra(w, f):
    nb,d = np.histogram(w,args.nbin)
    wb,d = np.histogram(w,args.nbin, weights=w)
    fb,d = np.histogram(w,args.nbin, weights=f)
    ok = nb > 0
    wb, fb = wb[ok] / nb[ok], fb[ok] / nb[ok] 
    fhi = 1.1 * np.percentile(fb, 99.9)
	
    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b')
    title = 'H-alpha lines of ' + spec[spec.rfind('/') + 1:] + ', ' + pos
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()
    
def total_flux(wavelength):
    _index = np.nonzero(w > wavelength)[0][0]
    _range = range(_index - 23, _index + 23)
    
    t_flux = sum([(w[i] - w[i-1]) * f[i] for i in _range])
    t_var = sum([(1/iv[i]) * (w_d[i] ** 2) for i in _range if iv[i] != 0])
    return t_flux, t_var
    
for spec in sdss:
    #print(spec)
    hdul = fits.open(spec)
    head  = hdul[0].header
    ra  = head['PLUG_RA']/15.
    dec = head['PLUG_DEC']
    pos = subs.d2hms(ra, sep = ' ', dp = 2) + ' ' + \
		  subs.d2hms(dec, sep = ' ', dp = 1, sign = True)
    #print('Position = ',pos)

    #print(hdul[3].data)
   
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
    
    #for i in range(0, w.size - 1):
    #   print(w[i] * 1.0)

    t.add_row([total_flux(6561.14892578)[0], total_flux(6561.14892578)[1], \
			   total_flux(4862.712890625)[0], total_flux(4862.712890625)[1], \
			   total_flux(4686.82128906)[0], total_flux(4686.82128906)[1]])
			   
    '''#Plot the entire spectra
    plot_spectra(w,f)
    
    # Plot H-alpha lines
    plot_spectra(w[ha_index - 23:ha_index + 23], \
                 f[ha_index - 23:ha_index + 23])
    
    # Plot H-beta lines
    plot_spectra(w[hb_index - 25:hb_index + 25], \
                 f[hb_index - 25:hb_index + 25])
                 
    # Plot HeII lines
    plot_spectra(w[he2_index - 25:he2_index + 25], \
                 f[he2_index - 25:he2_index + 25])'''
print(t)
