#!/usr/bin/env python

usage = \
"""
Plots SDSS DR12 spectra
"""

import os, sys, argparse
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
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

for spec in sdss:
    hdul = fits.open(spec)
    head  = hdul[0].header
    ra  = head['PLUG_RA']/15.
    dec = head['PLUG_DEC']
    pos = subs.d2hms(ra, sep = ' ', dp = 2) + ' ' + \
		  subs.d2hms(dec, sep = ' ', dp = 1, sign = True)
    print('Position = ',pos)
	
	# Store the range of flux and wavelength in arrays
    #table_bLines = hdul[3].data  
    table = hdul[1].data
    w  = subs.vac2air(10. ** table['loglam'])
    f  = table['flux']
    iv = table['ivar']
    ok = iv > 0.
    w  = w[ok]
    f  = f[ok]
    
    # Find the H-apha in the array of wavelengths
    for i in range(w.size // 2, w.size - 1):
        if w[i] > 6562.66:
            h_alpha_pos = i
            break
            
    # Calculate the total flux of H-alpha
    total_flux = 0
    for i in range(h_alpha_pos - 23, h_alpha_pos + 23):
        total_flux += (w[i] - w[i - 1]) * f[i]
    print('Total flux of H-alpha = ', total_flux)
    
	# Plot the entire spectra
    nb,d = np.histogram(w,args.nbin)
    wb,d = np.histogram(w,args.nbin, weights = w)
    fb,d = np.histogram(w,args.nbin, weights = f)
    ok = nb > 0
    wb = wb[ok] / nb[ok]
    fb = fb[ok] / nb[ok]
    fhi = 1.1 * np.percentile(fb, 99.9)

    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b')
    title = spec[spec.rfind('/') + 1:] + ', ' + pos
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()
    
    # Plot of H-alpha lines
    w_p = w[h_alpha_pos - 23:h_alpha_pos + 23]
    f_p = f[h_alpha_pos - 23:h_alpha_pos + 23]
    nb,d = np.histogram(w_p,args.nbin)
    wb,d = np.histogram(w_p,args.nbin, weights=w_p)
    fb,d = np.histogram(w_p,args.nbin, weights=f_p)
    ok = nb > 0
    wb = wb[ok] / nb[ok]
    fb = fb[ok] / nb[ok]
    fhi = 1.1 * np.percentile(fb, 99.9)

    fig=plt.figure(figsize = (12, 6))
    plt.plot(wb, fb, 'b')
    title = spec[spec.rfind('/') + 1:] + ', ' + pos
    plt.title(title.replace('_', '-'))
    plt.xlabel('Wavelength [\AA]')
    plt.show()

