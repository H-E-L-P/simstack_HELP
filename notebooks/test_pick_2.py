import pdb
import numpy as np
import pandas as pd
import os
from utils import clean_args
import astropy.units as u
try:
    from simstack import PickledStacksReader, measure_cib
except:
    from simstack.simstack import PickledStacksReader, measure_cib

path_pickles = os.environ['PICKLESPATH']
path_maps    = os.environ['MAPSPATH']
path_catalogs= os.environ['CATSPATH']

#Location of the stacked parameter file

os.system('ls -d /Users/Steven/Documents/python_simstack/output/simstack_fluxes/test2_IRAC* > all_run.txt')
runs = np.loadtxt('all_run.txt', usecols=[0], dtype = str)

file_config = 'example.cfg'


gal = np.zeros([45,45])
S250 = np.zeros([45,45])
S350 = np.zeros([45,45])
S500 = np.zeros([45,45])
name_z = np.array([])

k = 0
for name in runs:
  n = 0
  print name+'/'
  print '---------'
  stacked_flux_densities = PickledStacksReader(name+'/',file_config)
  flux =  stacked_flux_densities.simstack_flux_array
  flux2 = np.squeeze(flux)
  print dir(stacked_flux_densities)
  flux250 = flux2[0,:]
  flux350 = flux2[1,:]
  flux500 = flux2[2,:]
  for i in range(0,np.size(flux250)):
      zn = stacked_flux_densities.z_nodes[i:i+2]
      z_suf = '{:.2f}'.format(zn[0])+'-'+'{:.2f}'.format(zn[1])
      arg = clean_args('z_'+z_suf+'__m_12p80_22p00_sf')
      ng = len(stacked_flux_densities.bin_ids[arg])
      if ng > 0:
        print ng, flux250[i], flux350[i],flux500[i], '{:.2f}'.format(zn[1])
        S250[n,k] = flux250[i]
        S350[n,k] = flux350[i]
        S500[n,k] = flux500[i]
        gal[n,k] = ng
        n = n + 1 
        if k == 0:
          name_z = np.append(name_z,'{:.2f}'.format(zn[1]))
  k = k + 1

name_z = name_z.astype(np.float) 

S250_m = np.mean(S250, axis=1)
S350_m = np.mean(S350, axis=1)
S500_m = np.mean(S500, axis=1)
gal_m = np.mean(gal, axis=1)

cor_250 = S250[np.size(S250_m)-1,:]
cor_350 = S350[np.size(S350_m)-1,:]
cor_500 = S500[np.size(S500_m)-1,:]



loc = '/Users/Steven/Documents/python_simstack/output/simstack_fluxes/cats/'
DAT = np.array([gal_m,S250_m,S350_m,S500_m,name_z]).T
header = 'gal  S250  S350  S500  name'
np.savetxt(loc+'layer_flux.cat', DAT,delimiter= " ", fmt ='%s',header = header, newline='\n')