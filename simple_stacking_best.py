import numpy as np
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

  return (data1d - model)/err1d
  #return np.abs(data1d - model)/err1d


name = 'p_map_250_beth_noise_2pix.fits'
loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'
hdu = fits.open(loc+name)
img = hdu[1].data
imap = np.ndarray.flatten(img['a'])
ierr = np.ndarray.flatten(img['b'])
imap = imap - np.mean(imap)

print np.mean(imap) 
print np.size(imap)
print np.mean(ierr)
print np.sqrt(np.sum(imap**2)/np.size(imap))

imap = imap - np.mean(imap) 
name = 'pm_250_beth_noise_2pix.fits'
loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'
hdu = fits.open(loc+name)
img = hdu[1].data
cfits_flat = np.ndarray.flatten(img['a'])

ngal = np.loadtxt(loc+'ns_beth.cat')
len_model = len(imap)

SOURCE_LIST = np.append(np.linspace(12.4,18,15),np.linspace(18.4,26.4,21))
run = np.size(SOURCE_LIST)-1

header = 'Xi'
for q in range(run):
	header = header + ' TF'+str(SOURCE_LIST[q+1])

TF = np.zeros([run,run+2])
chi = np.zeros(run)+1e9

j_ont = 0   
#for k in range(2,run+2):
for k in range(30,run+2):

	t2 = time.time()
	cfits_flat_use = cfits_flat[:(k+1)*len_model]
 	nlayers = len(cfits_flat_use)/len_model 
 	chi_best = 1e9

	fit_params = Parameters()
	for iarg in range(nlayers):
		arg = 'name'+str(iarg)+'good'
		BG_use = -0.00758471515129
		if iarg == 0:
			#fit_params.add(arg,value= BG_use, min=BG_use-1e-12, max = BG_use+1e-12)  #, min=0.0, max = 1e-12)#, min=0.0
			fit_params.add(arg,value= 1e-3*np.random.randn())
		else:
			fit_params.add(arg,value= 1e-3*np.random.randn())


	cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
		  args=(cfits_flat_use,), kws={'data1d':imap,'err1d':ierr}, nan_policy = 'propagate') #, xtol = 1e-9, ftol = 1e-9

	values = np.array(cov_ss_1d.params)
	TF[k-2,0:k+1] = values
	model = np.zeros(len_model)
	for i in range(nlayers):
		model[:] += cfits_flat_use[i*len_model:(i+1)*len_model] * values[i]
	
	chi[k-2] = np.sum((imap - model)**2/ierr**2)
	LL = np.exp(-0.5*(imap - model)**2/ierr**2)/(np.sqrt(2*3.14159*ierr**2))
	bad = (imap - model)**2/ierr**2 < 6
	print k-2, chi[k-2], TF[k-2,0], np.prod(LL[bad])
	print LL[bad]
	print np.min(LL[bad]), np.max(LL[bad])
	exit()


 	t1 = time.time()
	tpass = t1-t2
	print "Total time run                       :",tpass/60, '  minuts'
exit()
DAT = (imap - model)/ierr
np.savetxt(loc+'test_r.cat', DAT.T, delimiter="   ", fmt='%s',  header= header, newline='\n')

exit()
DAT = np.zeros([run+2,run+1])
DAT[2:,0] = chi
DAT[:,1:run+1] = TF.T 
np.savetxt(loc+'res_beth_250_noise_2pix.cat', DAT, delimiter="   ", fmt='%s',  header= header, newline='\n')

t1 = time.time()
tpass = t1-t0
print "Total time                       :",tpass/60, '  minuts'

