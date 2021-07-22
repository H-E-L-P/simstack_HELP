#import pdb
import numpy as np
#import gc
#import os
#import os.path
#import sys
from utils import clean_args
from utils import clean_nans
from lmfit import Parameters, minimize #, fit_report
from astropy.io import fits
import time


t0 = time.time()


def simultaneous_stack_array_oned(p, layers_1d, data1d, err1d = None, arg_order = None):
  ''' Function to Minimize written specifically for lmfit '''
  v = p.valuesdict()
  len_model = len(data1d)
  nlayers = len(layers_1d)/len_model
  model = np.zeros(len_model)
  for i in range(nlayers):
  	model[:] += layers_1d[i*len_model:(i+1)*len_model] * v[v.keys()[i]]

  #if err1d is None:
  #return (data1d - model)/(np.sqrt(data1d)) 
  return (data1d - model)/err1d


#name = 'p_map_250_beth_noise.fits'
name = 'p_map_250_beth_noise_2pix.fits'
loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'
hdu = fits.open(loc+name)
img = hdu[1].data
imap = np.ndarray.flatten(img['a'])
ierr = np.ndarray.flatten(img['b'])
imap = imap - np.mean(imap)

#name = 'pm_250_beth.fits'
name = 'pm_250_beth_noise_2pix.fits'
loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'
hdu = fits.open(loc+name)
img = hdu[1].data
cfits_flat = np.ndarray.flatten(img['a'])
#aaa = cfits_flat < 0
#cfits_flat[aaa] = 0

ngal = np.loadtxt(loc+'ns_beth.cat')
len_model = len(imap)

SOURCE_LIST = np.append(np.linspace(12.4,18,15),np.linspace(18.4,26.4,21))
run = np.size(SOURCE_LIST)-1

header = 'BG '
for q in range(run):
	header = header + ' chi'+str(SOURCE_LIST[q+1])

for q in range(run):
	header = header + ' TF'+str(SOURCE_LIST[q+1])

#run_chi = 51
run_chi = 241
TF = np.zeros([run_chi,run])
chi = np.zeros([run_chi,run])+1e9
BG_use = np.linspace(-0.10,0.02,run_chi)
print BG_use

j_ont = 0   
mrun = run_chi+0
for k in range(2,run+2):
	t2 = time.time()
	cfits_flat_use = cfits_flat[:(k+1)*len_model]
 	nlayers = len(cfits_flat_use)/len_model 
 	chi_best = 1e9

	for j in range(j_ont,mrun):
		fit_params = Parameters()
		for iarg in range(nlayers):
			arg = 'name'+str(iarg)+'good'
			if iarg == 0:
				fit_params.add(arg,value= BG_use[j], min=BG_use[j]-1e-12, max = BG_use[j]+1e-12)  #, min=0.0, max = 1e-12)#, min=0.0
			else:
				fit_params.add(arg,value= 1e-3*np.random.randn())  #, min=0.0, max = 1e-12)#, min=0.0
		

		cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
		      args=(cfits_flat_use,), kws={'data1d':imap,'err1d':ierr}, nan_policy = 'propagate')

		values = np.array(cov_ss_1d.params)
		TF[j,k-2] = np.sum(ngal[:np.size(values)]*values)
		BG = values[0]
		model = np.zeros(len_model)
		for i in range(nlayers):
		  model[:] += cfits_flat_use[i*len_model:(i+1)*len_model] * values[i]

		#chi[j,k-2] = np.sum((imap - model)**2/imap)
		chi[j,k-2] = np.sum((imap - model)**2/ierr**2)
		print j, k-2, BG, chi[j,k-2], TF[j,k-2]
		if chi[j,k-2] < chi_best:
		  chi_best = chi[j,k-2]
		  j_ont = j

	mrun = np.min([run_chi,j_ont+30])
	j_ont = np.max([0,j_ont-30])
	print j_ont, mrun


 	t1 = time.time()
	tpass = t1-t2
	print "Total time run                       :",tpass/60, '  minuts'

	DAT = np.zeros([run_chi,run+run+1])
	DAT[:,0] = BG_use
	DAT[:,1:run+1] = chi 
	DAT[:,run+1:] = TF 
	#np.savetxt(loc+'res_beth_250_2.cat', DAT, delimiter="   ", fmt='%s',  header= header, newline='\n')
	np.savetxt(loc+'res_beth_250_noise_mean_2pixt.cat', DAT, delimiter="   ", fmt='%s',  header= header, newline='\n')

t1 = time.time()
tpass = t1-t0
print "Total time                       :",tpass/60, '  minuts'

