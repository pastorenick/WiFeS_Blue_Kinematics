#!/usr/bin/env python
# Filename: fit_vel_Nicola.py
# Run pPXF and MC simulations on single WiFeS spectra (blue arm)
#
# Output: dic_pPXF_results.dat (dictionary with spaxel coordinates, pPXF object, 
#         list of pPXF MC objects
#

#######################
# v.1 - it works!
# 
#######################


import pyfits, glob, numpy, os, time, pickle
from scipy import ndimage
from pylab import *
from astropy.io import ascii
from astropy.table import Table

from ppxf import ppxf
import ppxf_util as util


#################
### FUNCTIONS ###
#################

def get_wstart(ref, wave_ref, wave_per_pixel):
    return wave_ref - ((ref-1)*wave_per_pixel)

def get_wavelength(start_wave, wave_per_pixel, size):
    return np.array([start_wave + i*wave_per_pixel for i in range(size)])


def findCoords_from_string(filename):   #Retrieves the spaxel coordinates from the filename
  #
  listGood_indices, coordinates = [], []
  for ii in numpy.arange(len(filename)):
    if filename[-ii] == '_':
      try:
        tmp = float(filename[-ii+1])
        try: 
          tmp2 = float(filename[-ii+2])
          coordinates.append(tmp2+10*tmp) 
        except:
          coordinates.append(tmp)
      except:
        dummy = 1
  if numpy.isreal(coordinates[0]) and numpy.isreal(coordinates[1]): 
    return coordinates
  else:
    print "Fail"
    return False



def randomArr(sigmaGauss_arr): #Takes a sigma array and create a gaussian random array
  from random import gauss
  Random_arr = []
  for ii in sigmaGauss_arr:
    Random_arr.append(gauss(0., ii))
  return numpy.array(Random_arr)

##############
#### MAIN ####
##############

t1 = time.time()
# Read a galaxy spectrum and define the wavelength range
#
directory = './'
filenames = glob.glob('./Blue_Ascii_specs/*')   #Reads the ascii files with all the blue spectra

name, rec, sig, chi2, er = [], [], [], [], []

if not(os.path.exists('./Plot')): os.mkdir('./Plot') 
     

dic_pPXF_results = {}

moments = 2
N_iter_MC = 100


for ii in numpy.arange(len(filenames)):
  print ii, "/", len(filenames)
  wave, gal_lin, variance_lin, pix = numpy.genfromtxt(filenames[ii], unpack = True, skiprows = 1) 
  Spaxel_coordinates = findCoords_from_string(filenames[ii])
  #
  text =  "Blue_Spec_"+str(int(Spaxel_coordinates[0]))+"_"+str(int(Spaxel_coordinates[1]))
  #
  noise_lin = numpy.sqrt(variance_lin)
  out_range = numpy.where(numpy.logical_or(wave< 4020., wave > 5550.))
  wave_sel = numpy.delete(wave, out_range)
  gal_lin_sel = numpy.delete(gal_lin, out_range)
  noise_lin_sel = numpy.delete(noise_lin, out_range)
  lamRange1 = numpy.array([wave_sel[0], wave_sel[-1]])
  #
  galaxy, logLam1, velscale = util.log_rebin(lamRange1, gal_lin_sel)
  fac = numpy.median(galaxy)
  galaxy = galaxy / fac# Normalize spectrum to avoid numerical issues
  error, logLam1, velscale = util.log_rebin(lamRange1, noise_lin_sel**2)
  noise = numpy.sqrt(error)
  noise = noise / fac
  #
  #
  logLam2, fl_temp = numpy.genfromtxt(directory + 'wise_template.dat', unpack = True, skiprows = 1)
  wave_temp = numpy.exp(logLam2)
  lamRange2 = numpy.array([logLam2[0], logLam2[-1]])
  #
  templates = fl_temp
  #
  c = 299792.458
  dv = (logLam2[0]-logLam1[0])*c # km/s
  #
  vel = 1500. 
  #
  start = [vel, 200.] # (km/s), starting guess for [V,sigma]
  #
  pp = ppxf(templates, galaxy, noise, velscale, start, plot=False, moments=2, 
        degree=4, vsyst=dv)
  #
  if moments == 4:
    pp = ppxf(templates, galaxy, noise, velscale, [pp.sol[0], pp.sol[1]], plot=False, moments=2, 
        degree=4, vsyst=dv)
  #
  #
  rec.append(pp.sol[0])
  sig.append(pp.sol[1])
  chi2.append(pp.chi2)
  er.append(pp.error*np.sqrt(pp.chi2))
  #
  # RUN MONTECARLO SIMULATION
  #
  lenMC = N_iter_MC
  #
  listPP_MC = []
  for jj in numpy.arange(lenMC):
    # Create fake spectrum
    fake_galaxy = pp.bestfit+randomArr(noise)
    # Run pPXF on fake spectrum
    pp_MC = ppxf(templates, fake_galaxy, noise, velscale, start, plot=False, 
        moments=2, degree=4, vsyst=dv, quiet=True)
    #
    if moments == 4:
      pp = ppxf(templates, galaxy, noise, velscale, [pp.sol[0], pp.sol[1]], plot=False, moments=2, 
        degree=4, vsyst=dv, quiet=True) 
    # Store pPXF fitting result
    listPP_MC.append(pp_MC.sol)
  #
  dic_pPXF_results[ii] = [Spaxel_coordinates, pp, listPP_MC, 'Coordinates, pp object, pp_MC.sol']
  #

fileOut = open('./dic_pPXF_results.dat', 'wb')
pickle.dump(dic_pPXF_results, fileOut)
fileOut.close()

print "#########################"
print "Done in "+str(round(time.time()-t1),2)+"s"
print "########################"
