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
 
  #model -= np.mean(model) ###removed by me!

  #return (data1d - model)/err1d
  return ( (data1d - model)/err1d )


loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'
#ngal = np.loadtxt(loc+'ns_beth.cat')
#ngal = np.loadtxt(loc+'ns_beth_2arc_30arc_4ac.cat')
#ngal = np.loadtxt(loc+'ns_beth_tiny01.cat')

#ngal[0:np.size(ngal)-13] =0


back = 0 # -1*np.mean(imap)

#should = np.loadtxt(loc+'should.cat')
#should = np.append([back,0],should)

#print np.sum(ngal*should)
#-------------
#res = np.loadtxt(loc+'test_r.cat')
#good = res**2 < 25
#nlayers = len(cfits_flat)/len_model 
#cfits_flat_new = np.array([])
#for iarg in range(nlayers):
#  val = cfits_flat[iarg*len_model:(iarg+1)*len_model]
#  aaa = val < 1e-9
#  print np.size(val[aaa])
#  print 'iarg',np.mean(val), np.max(val), np.std(val), np.min(val)

  #use = cfits_flat[iarg*len_model:(iarg+1)*len_model]
  #print np.size(use), np.size(imap)
  #cfits_flat_new = np.append(cfits_flat_new,use[good])

#cfits_flat = cfits_flat_new
#imap = imap[good]
#imap = imap - np.mean(imap)
#ierr = ierr[good]
#len_model = len(imap)
#-------------


SOURCE_LIST = np.append(np.linspace(12.4,18,15),np.linspace(18.4,26.4,21))
'''
SOURCE_LIST = np.linspace(12.5,26.5,15)
SOURCE_LIST = np.linspace(12.5,26.5,29)
SOURCE_LIST = np.linspace(12.5,26.5,8)
SOURCE_LIST = np.linspace(12.4,26.4,71)
SOURCE_LIST = np.linspace(12.4,26.4,2)
SOURCE_LIST = np.linspace(12.4,26.4,141)
'''
run = np.size(SOURCE_LIST)-1

header = 'Xi'
for q in range(run):
	header = header + ' TF'+str(SOURCE_LIST[q+1])

boot = 10
TF = np.zeros([boot,run+2])
chi = np.zeros(boot)+1e9
total = np.zeros(10)

#names = ['full','full_noise','full_BG','full_noise_BG','full_mean','full_mean_noise','full_mean_BG','full_mean_noise_BG', \
names = ['full_mean_noise_BG','pois_full','pois_full_noise','pois_full_BG','pois_full_noise_BG','pois_full_mean','pois_full_mean_noise','pois_full_mean_BG','pois_full_mean_noise_BG']

names = ['full_tiny01_BG_001','full_tiny01','full_tiny01_BG','full_pois_tiny01','full_pois_tiny01_BG']
names = ['full_fwhm1','full_fwhm1_BG']
names = ['full_BG_276']
names = ['full_232plus05mjy_BG_276','full_232plus01mjyA_BG_276']
names = ['full_232plus03mjyA_BG_276']
names = ['full_BG_276_rd']
names = ['full_noise_2arc_30arc_BG_4ac']
names = ['full_noise_simstack']

names = ['r','K','IRAC']#,'MIPS','PACS','SPIRE','S2CLS','VLA']#,'EN1','CDFS']

redshift_nodes  = np.array([0.01, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,10]) 


#names = ['all_1.77_mag_MY0.5','all_1.77_mag_MY1.0','all_1.77_mag_MY1.5','all_1.77_mag_MY2.0','all_1.77_mag_MY2.5'\
#		,'all_1.77_mag_MY3.0','all_1.77_mag_MY3.5', 'all_1.77_mag_MY4.0','all_1.77_mag_MY10.0'\
names =['IRAC2_MY0.5','IRAC2_MY1.0','IRAC2_MY1.5','IRAC2_MY2.0','IRAC2_MY2.5'\
		,'IRAC2_MY3.0','IRAC2_MY3.5', 'IRAC2_MY4.0','IRAC2_MY10.0'\
        ,'r2_MY0.5','r2_MY1.0','r2_MY1.5','r2_MY2.0','r2_MY2.5'\
		,'r2_MY3.0','r2_MY3.5', 'r2_MY4.0','r2_MY10.0','K2_MY0.5','K2_MY1.0','K2_MY1.5','K2_MY2.0','K2_MY2.5'\
		,'K2_MY3.0','K2_MY3.5', 'K2_MY4.0','K2_MY10.0']

#names = ['r2_4ac','K2_4ac','IRAC2_4ac','MIPS_4ac','PACS_4ac','SPIRE_4ac','S2CLS_4ac','VLA_4ac']	
#names = ['IRAC2_6ac','IRAC2_2ac','IRAC2_4ac']	
#names = ['all_1.38_mag_MY0.5','all_1.38_mag_MY1.0','all_1.38_mag_MY1.5','all_1.38_mag_MY2.0','all_1.38_mag_MY2.5'\
#		,'all_1.38_mag_MY3.0','all_1.38_mag_MY3.5', 'all_1.38_mag_MY4.0','all_1.38_mag_MY10.0']

names = ['MIPS','PACS','SPIRE','S2CLS','VLA']	

#names = ['IRAC2_cv','EN1_cv','CDFS_cv']
names = ['all_1.38_mag_MY0.5','all_1.38_mag_MY1.0','all_1.38_mag_MY1.5','all_1.38_mag_MY2.0',\
		'all_1.38_mag_MY2.5','all_1.38_mag_MY3.0','all_1.38_mag_MY3.5',\
         'all_1.38_mag_MY4.0','all_1.38_mag_MY10.0']
names = ['all_1.38_mag_4ac','all_1.38_mag'] 
names = ['all_1.38_mag_4ac']        

names = ['full_noise_2arc_30arc2ac_MY0.5','full_noise_2arc_30arc2ac_MY1.0','full_noise_2arc_30arc2ac_MY1.5',\
		'full_noise_2arc_30arc2ac_MY2.0','full_noise_2arc_30arc2ac_MY2.5','full_noise_2arc_30arc2ac_MY3.0',\
		'full_noise_2arc_30arc2ac_MY3.5','full_noise_2arc_30arc2ac_MY4.0','full_noise_2arc_30arc2ac_MY10.0']#,'full_noise_simstack3.5']
#names = ['r2','K2','IRAC2']
names = ['all_1.38_mag_5','all_1.38_mag_5_4ac'] 

names = ['r2_FWHM5.0']
names = ['K_uds24']
names = ['H_huble_JK1','H_huble_JK2','H_huble_JK3','H_huble_JK4']
#names = ['K_uds_JK1','K_uds_JK2','K_uds_JK3','K_uds_JK4']

#names = ['r2_JK4','K2_JK4','IRAC2_JK4','all_1.38_mag_3_JK4']
#names = ['MIPS_JK4','PACS_JK4','SPIRE_JK4','S2CLS_JK4','VLA_JK4']	

names = ['r2_PSF20.0','r2_PSF30.0','r2_PSF20.0_4ac','r2_PSF30.0_4ac','r2_PSF40.0','r2_PSF40.0_4ac']
names = ['r2_PSF19.0','r2_PSF19.0_4ac']
names = ['r2_FWHM17.0','r2_FWHM17.0_4ac']

wave = '250'
for q in range(0,np.size(names)):
#for q in range(3,4):

	#ngal = np.loadtxt(loc+'ns_beth_2arc_30arc.cat')

	name2 = names[q]
	print name2
	#name = 'p_map_'+wave+'_beth_'+name2+'.fits'

	#name = 'p_map_'+wave+'_'+name2+'_MY.fits'
	name = 'p_map_'+wave+'_'+name2+'.fits'
	hdu = fits.open(loc+name)
	img = hdu[1].data
	imap = np.ndarray.flatten(img['a'])
	ierr = np.ndarray.flatten(img['b'])


	imap = imap - np.mean(imap)
	#name = 'pm_'+wave+'_beth_'+name2+'.fits'

	#name = 'pm_'+wave+'_beth_'+name2+'_MY.fits'
	name = 'pm_'+wave+'_beth_'+name2+'.fits'
	hdu = fits.open(loc+name)
	img = hdu[1].data
	cfits_flat = np.ndarray.flatten(img['a'])


	len_model = len(imap)

	print 'helo'
	t2 = time.time()
	nlayers = len(cfits_flat)/len_model 
	chi_best = 1e9

	cfits_flat_use = np.zeros(np.size(cfits_flat))

	bs_sample = np.random.randint(0, len_model,len_model)

	fit_params = Parameters()
	for iarg in range(nlayers):
	  arg = 'name'+str(iarg)+'good'
	  fit_params.add(arg,value= 1e-6*np.random.randn())
	  ## coment out
	  use = cfits_flat[iarg*len_model:(iarg+1)*len_model]
	  cfits_flat_use[iarg*len_model:(iarg+1)*len_model] = use[bs_sample]

	#cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
	#	  args=(cfits_flat_use,), kws={'data1d':imap[bs_sample],'err1d':ierr[bs_sample]}, nan_policy = 'propagate')

	cov_ss_1d = minimize(simultaneous_stack_array_oned, fit_params,
		  args=(cfits_flat,), kws={'data1d':imap,'err1d':ierr}, nan_policy = 'propagate')#, xtol = 1e-15, ftol = 1e-15, epsfcn = 0.001) 


	values = np.array(cov_ss_1d.params)
	#TF[q,:] = values
	print(values)
	model = np.zeros(len_model)
	for i in range(nlayers):
		model[:] += cfits_flat[i*len_model:(i+1)*len_model] * values[i]
	
	#chi[q] = np.sum((imap - model)**2/ierr**2)
	#print q, chi[q], values[0], np.sum(ngal*values),cov_ss_1d.redchi, cov_ss_1d.chisqr, cov_ss_1d.success, cov_ss_1d.nfev

	##model = np.zeros(len_model)
	##for i in range(nlayers):
	##	model[:] += cfits_flat[i*len_model:(i+1)*len_model] * should[i]


 	t1 = time.time()
	tpass = t1-t2	
	print "Total time run                       :",tpass/60, '  minuts'
    
	DAT =  np.array(values).T
	##np.savetxt(loc+'res_beth_'+wave+'_'+name2+'.cat', DAT, delimiter="   ", fmt='%s', newline='\n')
	np.savetxt(loc+'res_'+wave+'_'+name2+'.cat', DAT, delimiter="   ", fmt='%s', newline='\n')

exit()

DAT = np.zeros([run+2,boot])
#DAT[:,0] = chi
DAT[:,:] = TF.T 
np.savetxt(loc+'res_beth_250_noise_2pix_BOOT_25.cat', DAT, delimiter="   ", fmt='%s',  header= header, newline='\n')

t1 = time.time()
tpass = t1-t0
print "Total time                       :",tpass/60, '  minuts'
#print np.mean(total)
#print np.std(total)


