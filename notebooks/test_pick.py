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
shortnmae = 'test'
path_config = path_pickles + '/output/simstack_fluxes/test/'
path_config = path_pickles + '/output/bootstrapped_fluxes/test/'

#file_config = 'uvista__DR2__2pop__7_maps_s15_binning.cfg'
#shortname = 'uVista_Laigle_v1.1__2pop__7bands__s15_bins_in_slices'
shortname = 'IRAC2_test_all_z'
shortname =  'MOCK_r_test_all_z_boot_0'
##shortname =  'MOCK_r_c_test_all_z'
##shortname =  'MIPS_test_all_z'
file_config = 'example2.cfg'


#print path_config
#print file_config
stacked_flux_densities = PickledStacksReader(path_config,file_config)


#print stacked_flux_densities.params.keys()
print np.shape(stacked_flux_densities.bootstrap_flux_array)
#print stacked_flux_densities.params['bins']


names = ['r','K','IRAC','MIPS','PACS','SPIRE','S2CLS','VLA','MOCKr','IRAC_extra_layer','IRAC_not_mean','IRAC_no_negative_layer', \
           'IRAC_no_negative','IRAC_no_negative_tmap','IRAC_tmap']
#names = ['rmock','r_mc_old','K_mc_old','IRAC_mc_old','MOCKr_extra_sources']           
names = names[8]


#print dir(stacked_flux_densities)
#print stacked_flux_densities.bin_ids
#print np.shape(stacked_flux_densities.simstack_nuInu_array)

##flux =  stacked_flux_densities.simstack_flux_array
flux =  stacked_flux_densities.bootstrap_flux_array
err = stacked_flux_densities.boot_error_bars

err2 = np.squeeze(err)
flux2 = np.squeeze(flux)

flux250 = flux2[0,:]
flux350 = flux2[1,:]
flux500 = flux2[2,:]

#print flux2[1,:]

gal = np.array([])
S250 = np.array([])
S350 = np.array([])
S500 = np.array([])
name = np.array([])

#print stacked_flux_densities.bin_ids
print 'begin'
#for i in range(stacked_flux_densities.nz):
for i in range(0,np.size(flux250)):
    print stacked_flux_densities.bin_ids
    exit()
    zn = stacked_flux_densities.z_nodes[i:i+2]
    z_suf = '{:.2f}'.format(zn[0])+'-'+'{:.2f}'.format(zn[1])
    arg = clean_args('z_'+z_suf+'__m_12p80_22p00_sf')
    ng = len(stacked_flux_densities.bin_ids[arg])
    print ng
    if ng > 0:
      print ng, flux250[i], flux350[i],flux500[i], '{:.2f}'.format(zn[1])
      S250 = np.append(S250,flux250[i])
      S350 = np.append(S350,flux350[i])
      S500 = np.append(S500,flux500[i])
      gal = np.append(gal,ng)
      name = np.append(name,'{:.2f}'.format(zn[1]))
name =  name.astype(np.float)

loc = '/Users/Steven/Documents/python_simstack/output/simstack_fluxes/cats/'
DAT = np.array([gal,S250,S350,S500,name]).T
header = 'gal  S250  S350  S500  name'
np.savetxt(loc+names+'_flux_boot.cat', DAT,delimiter= " ", fmt ='%s',header = header, newline='\n')




