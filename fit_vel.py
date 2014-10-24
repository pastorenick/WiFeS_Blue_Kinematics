#!/usr/bin/env python
##############################################################################
#
# Usage example for the procedure PPXF, which
# implements the Penalized Pixel-Fitting (pPXF) method by
# Cappellari M., & Emsellem E., 2004, PASP, 116, 138.
# The example also shows how to include a library of templates
# and how to mask gas emission lines if present.
# 
# MODIFICATION HISTORY:
#   V1.0: Written by Michele Cappellari, Leiden 11 November 2003
#   V1.1: Log rebin the galaxy spectrum. Show how to correct the velocity
#       for the difference in starting wavelength of galaxy and templates.
#       MC, Vicenza, 28 December 2004
#   V1.11: Included explanation of correction for instrumental resolution.
#       After feedback from David Valls-Gabaud. MC, Venezia, 27 June 2005
#   V2.0: Included example routine to determine the goodPixels vector
#       by masking known gas emission lines. MC, Oxford, 30 October 2008
#   V2.01: Included instructions for high-redshift usage. Thanks to Paul Westoby
#       for useful feedback on this issue. MC, Oxford, 27 November 2008
#   V2.02: Included example for obtaining the best-fitting redshift.
#       MC, Oxford, 14 April 2009
#   V2.1: Bug fix: Force PSF_GAUSSIAN to produce a Gaussian with an odd 
#       number of elements centered on the middle one. Many thanks to 
#       Harald Kuntschner, Eric Emsellem, Anne-Marie Weijmans and 
#       Richard McDermid for reporting problems with small offsets 
#       in systemic velocity. MC, Oxford, 15 February 2010
#   V2.11: Added normalization of galaxy spectrum to avoid numerical
#       instabilities. After feedback from Andrea Cardullo.
#       MC, Oxford, 17 March 2010
#   V2.2: Perform templates convolution in linear wavelength. 
#       This is useful for spectra with large wavelength range. 
#       MC, Oxford, 25 March 2010
#   V2.21: Updated for Coyote Graphics. MC, Oxford, 11 October 2011
#   V2.22: Renamed PPXF_KINEMATICS_EXAMPLE_SAURON to avoid conflict with the 
#       new PPXF_KINEMATICS_EXAMPLE_SDSS. Removed DETERMINE_GOOPIXELS which was 
#       made a separate routine. MC, Oxford, 12 January 2012
#   V3.0: Translated from IDL into Python. MC, Oxford, 6 December 2013
#       
##############################################################################
import pyfits
from scipy import ndimage
import numpy as np
from time import clock
import glob
import pylab as pl
from astropy.io import ascii
from astropy.table import Table
import os


from ppxf import ppxf
import ppxf_util as util

def get_wstart(ref, wave_ref, wave_per_pixel):
    return wave_ref - ((ref-1)*wave_per_pixel)

def get_wavelength(start_wave, wave_per_pixel, size):
    return np.array([start_wave + i*wave_per_pixel for i in range(size)])

def ppxf_kinematics_example_sauron():

    # Read a galaxy spectrum and define the wavelength range
    #
    dir = 'temp_wise/'
    #file = '/Users/martinaf/Documents/work/Swinburne/WISE/stacked_spectrumBlue.dat'
    file = glob.glob('/Users/martinaf/Documents/work/Swinburne/WISE/Blue_Ascii_specs/*')

    name = []
    rec = []
    sig = []
    chi2 = []
    er = []


    if not(os.path.exists('../WISE/Plot')):
        os.mkdir('../WISE/Plot')

    for ii in range(len(file)):
        print ii
        wave, gal_lin, noise_lin, pix = np.genfromtxt(file[ii], unpack = True, skiprows = 1) 
        text =  file[ii][63:80]
        name.append(text)

        noise_lin = np.sqrt(noise_lin)
        out_range = np.where(np.logical_or(wave< 4020., wave > 5550.))
        wave = np.delete(wave, out_range)
        gal_lin = np.delete(gal_lin, out_range)
        noise_lin = np.delete(noise_lin, out_range)
        lamRange1 = np.array([wave[0], wave[-1]])

        galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal_lin)
        fac = np.median(galaxy)
        galaxy = galaxy / fac# Normalize spectrum to avoid numerical issues
        error, logLam1, velscale = util.log_rebin(lamRange1, noise_lin**2)
        noise = np.sqrt(error)
        noise = noise / fac


        logLam2, fl_temp = np.genfromtxt(dir + 'wise_template.dat', unpack = True, skiprows = 1)
        wave_temp = np.exp(logLam2)
        lamRange2 = np.array([logLam2[0], logLam2[-1]])

        templates = fl_temp

        c = 299792.458
        dv = (logLam2[0]-logLam1[0])*c # km/s
    
        vel = 1500. 

        start = [vel, 200.] # (km/s), starting guess for [V,sigma]
        t = clock()

        pp = ppxf(templates, galaxy, noise, velscale, start, plot=False, moments=2, 
              degree=4, vsyst=dv)

        rec.append(pp.sol[0])
        sig.append(pp.sol[1])
        chi2.append(pp.chi2)
        er.append(pp.error*np.sqrt(pp.chi2))

        # NICOLA
        # fig = pl.figure(0)
        # pl.clf()
        # ax1 = pl.subplot(111)
        # ax1.plot(np.exp(logLam1), pp.galaxy, 'k-')
        # ax1.plot(np.exp(logLam1), pp.bestfit, 'r-')
        # pl.savefig('../WISE/Plot/'+text[:-3]+'.pdf', bbox_inches='tight')
        # END NICOLA



    data = Table([name, rec, sig, chi2, er], names=['position', 'velocity', 'sigma', 'chi2', 'error'])
    ascii.write(data, '/Users/martinaf/Documents/work/Swinburne/WISE/velocity1')

        # par = open('/Users/martinaf/Documents/work/Swinburne/WISE/velocity'+text,'w')
        # par.write(str(pp.sol[0])+ ' ' + str(pp.sol[1]))
        # par.close()

    
    # print "Formal errors:"    
    # print "     dV    dsigma   dh3      dh4"
    # print "".join("%8.2g" % f for f in pp.error*np.sqrt(pp.chi2))

    # print 'elapsed time in PPXF (s):', clock() - t 

    # logLam1 = np.exp(logLam1)

    # resid = pp.galaxy - pp.bestfit

    # nPixels = pp.galaxy.shape[0] #va cambiato se lo fai per il bestfit

    # x = np.linspace(-1, 1, len(pp.galaxy))
    # apoly = np.polynomial.legendre.legval(x, pp.polyweights)


    # red = pp.sol[0]/c

    # logLam1 = logLam1 / (1 + red)
    # pp.galaxy = pp.galaxy * (1 + red)
    # pp.bestfit = pp.bestfit * (1 + red)
    # pp.noise = pp.noise * (1 + red)
    # resid = resid * (1 + red)
    # apoly = apoly * (1 + red)




         
#------------------------------------------------------------------------------

if __name__ == '__main__':
    ppxf_kinematics_example_sauron()
