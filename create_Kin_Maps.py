#!/usr/bin/env python
# Filename: create_Kin_Maps.py
# Reads the output from fit_vel_Nicola.py and creates the velocity and 
# sigma maps for N1407
#
# Output: 

#######################
# v.0.1 - 
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



def drawMap(arr2d, interac = True, zunity = 'counts', axx = 'None', zlimits = None):  #Takes a 2d array (e.g. 38x25 arcsec) and plot it
  if axx == 'None':
    if interac: plt.ion()
    fig = figure(num=0)
    axmap1 = subplot(111)
    if zlimits == None:
      map1 = axmap1.imshow(arr2d, origin='lower', interpolation='nearest', aspect='equal', 
                vmin = numpy.min(arr2d), vmax=numpy.max(arr2d))
    else:
      map1 = axmap1.imshow(arr2d, origin='lower', interpolation='nearest', aspect='equal', 
                vmin = zlimits[0], vmax=zlimits[1])
    cbar = fig.colorbar(map1, ax = axmap1, orientation='vertical', shrink=1)
    axmap1.set_xlabel('X [arcsec]')
    axmap1.set_ylabel('Y [arcsec]')
    cbar.set_label('Z ['+zunity+']')
    return fig
  else:
    if zlimits == None:
      map1 = axx.imshow(arr2d, origin='lower', interpolation='nearest', aspect='equal', 
                vmin = numpy.min(arr2d), vmax=numpy.max(arr2d))
    else:
      map1 = axx.imshow(arr2d, origin='lower', interpolation='nearest', aspect='equal', 
                vmin = zlimits[0], vmax=zlimits[1])
    colorbar(map1, ax = axx, orientation='vertical', shrink=1)
    axx.set_xlabel('X [arcsec]')
    axx.set_ylabel('Y [arcsec]')
    return True

##############
#### MAIN ####
##############

fileIn = open('./dic_pPXF_results.dat', 'wb')
dic_pPXF_results = pickle.load(fileOut)
fileIn.close()


[Spaxel_coordinates, pp, listPP_MC, 'Coordinates, pp object, pp_MC object'] = 

Moments = len(dic_pPXF_results[1].sol)






  header['RA'] = '03:40:11.860' 
  header['Dec'] = '-18:34:48.4'


fileIn = open('pPXFresult4.dat', 'rb')
solutionsPpxf = pickle.load(fileIn)
fileIn.close()

#Reshaping map
mapKin_Bin_plain = numpy.empty((len(X.ravel()), 8))
mapKin_Bin = numpy.empty((len(X), len(X[0]), 8))


for kk in range(4):
  for ii in range(len(elementsVor)):
    for jj in elementsVor[ii]:
      mapKin_Bin_plain[jj,kk] = solutionsPpxf[ii][kk] #Vel disp
  mapKin_Bin[:,:, kk] = mapKin_Bin_plain[:,kk].reshape(shape(X))

fig = figure(figsize=(10,12))

ax = []

zlim = [[-45,45], [120,300], [-0.1,0.1], [-0.1,0.1]]  #Following Jacob's plots

for kk in range(4):
  ax.append(subplot(2,2,kk+1))
  drawMap(mapKin_Bin[:, :, kk], zunity = 'km/s', axx = ax[kk], zlimits = zlim[kk])







fig = figure(num=5, figsize=(12,10))
ax1 = subplot(221)
ax1.plot(elldistBin[indices], velBin[indices], 'ko')
ax1.set_title('Velocity')
ax1.set_xlabel('Galactocentric distance [arcsec]')
ax1.set_xlabel('Velocity [km/s]')

ax2 = subplot(222)
ax2.plot(elldistBin[indices], sigmaBin[indices], 'ko')
ax2.set_title('Velocity Dispersion')
ax2.set_xlabel('Galactocentric distance [arcsec]')
ax2.set_xlabel('Velocity dispersion [km/s]')


ax3 = subplot(223)
ax3.plot(elldistBin[indices], h3Bin[indices], 'ko')
ax3.set_title('h3')
ax3.set_xlabel('Galactocentric distance [arcsec]')
ax3.set_xlabel('h3')


ax4 = subplot(224)
ax4.plot(elldistBin[indices], h4Bin[indices], 'ko')
ax4.set_title('h4')
ax4.set_xlabel('Galactocentric distance [arcsec]')
ax4.set_xlabel('h4 [km/s]')

savefig('Results/Kinematics.pdf', bbox_inches='tight')
