from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from utils import gauss_kern
from utils import circle_mask
from utils import clean_nans
from utils import smooth_psf
from astropy.table import Table
import scipy.signal

from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as units
from astropy.coordinates import match_coordinates_sky

from sklearn.utils import shuffle
import pandas as pd

loc_m = '/Users/Steven/Documents/python_simstack/maps/'

################################################### names of files
list_names = ['MIPS','PACS','SPIRE','S2CLS','VLA']#,'EN1','CDFS']'r2','K2','IRAC2',
#map_names = ['cosmos_4.0_arcsec_pixels_PSW.fits',\
#             'cosmos_4.0_arcsec_pixels_PSW.fits','cosmos_4.0_arcsec_pixels_PSW.fits','cosmos_4.0_arcsec_pixels_PSW.fits',\
#             'cosmos_4.0_arcsec_pixels_PSW.fits']#,'SPIRE250_mask_EN1_croped.fits','SPIRE250_err_CDFS.fits'] 'SPIRE250_mask4.fits','SPIRE250_mask4K.fits','SPIRE250_mask4.fits',

#map_names = ['cosmos_4.0_arcsec_pixels_PMW.fits',\
#             'cosmos_4.0_arcsec_pixels_PMW.fits','cosmos_4.0_arcsec_pixels_PMW.fits','cosmos_4.0_arcsec_pixels_PMW.fits',\
#             'cosmos_4.0_arcsec_pixels_PMW.fits']#,'SPIRE350_mask_EN1_croped.fits','SPIRE350_err_CDFS.fits'] SPIRE350_mask4.fits','SPIRE350_mask4K.fits','SPIRE350_mask4.fits'

map_names = ['cosmos_4.0_arcsec_pixels_PLW.fits',\
             'cosmos_4.0_arcsec_pixels_PLW.fits','cosmos_4.0_arcsec_pixels_PLW.fits','cosmos_4.0_arcsec_pixels_PLW.fits',\
             'cosmos_4.0_arcsec_pixels_PLW.fits']#,'SPIRE500_mask_EN1_croped.fits','SPIRE500_err_CDFS.fits']''' 'SPIRE500_mask4.fits','SPIRE500_mask4K.fits','SPIRE500_mask4.fits',

#emap_names = ['cosmos_4.0_arcsec_pixels_PSW.fits',\
#             'cosmos_4.0_arcsec_pixels_PSW.fits','cosmos_4.0_arcsec_pixels_PSW.fits','cosmos_4.0_arcsec_pixels_PSW.fits',\
#             'cosmos_4.0_arcsec_pixels_PSW.fits'] #'SPIRE250_err4.fits','SPIRE250_err4.fits','SPIRE250_err4.fits',

#emap_names = ['cosmos_4.0_arcsec_pixels_PMW.fits',\
#             'cosmos_4.0_arcsec_pixels_PMW.fits','cosmos_4.0_arcsec_pixels_PMW.fits','cosmos_4.0_arcsec_pixels_PMW.fits',\
#             'cosmos_4.0_arcsec_pixels_PMW.fits'] #'SPIRE350_err4.fits','SPIRE350_err4.fits','SPIRE350_err4.fits',

emap_names = ['cosmos_4.0_arcsec_pixels_PLW.fits',\
             'cosmos_4.0_arcsec_pixels_PLW.fits','cosmos_4.0_arcsec_pixels_PLW.fits','cosmos_4.0_arcsec_pixels_PLW.fits',\
             'cosmos_4.0_arcsec_pixels_PLW.fits']  #'SPIRE500_err4.fits','SPIRE500_err4.fits','SPIRE500_err4.fits',


'''
list_names = ['r2','K2','IRAC2','all_1.77_mag']
map_names = ['SPIRE250_mask4.fits','SPIRE250_mask4K.fits','SPIRE250_maskK.fits','SPIRE250_mask4.fits']
emap_names = ['SPIRE250_err4.fits','SPIRE250_err4.fits','SPIRE250_err4.fits','SPIRE250_err4.fits']

list_names = ['IRAC2','all_1.77_mag']
map_names = ['SPIRE250_mask4.fits','SPIRE250_mask4K.fits']
emap_names = ['SPIRE250_err4.fits','SPIRE250_err4.fits']


list_names = ['r2','K2','IRAC2','all_1.38_mag_3']#,'all_1.77_mag']
map_names = ['SPIRE250_mask4.fits','SPIRE250_mask4K.fits','SPIRE250_mask4.fits','SPIRE250_mask4K.fits']
emap_names = ['SPIRE250_err4.fits','SPIRE250_err4.fits','SPIRE250_err4.fits','SPIRE250_err4.fits']

map_names = ['SPIRE350_mask4.fits','SPIRE350_mask4K.fits','SPIRE350_mask4.fits','SPIRE350_mask4K.fits']
emap_names = ['SPIRE350_err4.fits','SPIRE350_err4.fits','SPIRE350_err4.fits','SPIRE350_err4.fits']

map_names = ['SPIRE500_mask4.fits','SPIRE500_mask4K.fits','SPIRE500_mask4.fits','SPIRE500_mask4K.fits']
emap_names = ['SPIRE500_err4.fits','SPIRE500_err4.fits','SPIRE500_err4.fits','SPIRE500_err4.fits']


list_names = ['IRAC2']
map_names = ['SPIRE500_mask4.fits']
emap_names = ['SPIRE500_err4.fits']

list_names = ['EN1','CDFS']
map_names = ['SPIRE250_mask_EN1_croped.fits','SPIRE250_mask_CDFS.fits']
emap_names = ['SPIRE250_err_EN1_croped.fits','SPIRE250_err_CDFS.fits']

#map_names = ['SPIRE350_mask_EN1_croped.fits','SPIRE350_mask_CDFS.fits']
#emap_names = ['SPIRE350_err_EN1_croped.fits','SPIRE350_err_CDFS.fits']

#map_names = ['SPIRE500_mask_EN1_croped.fits','SPIRE500_mask_CDFS.fits']
#emap_names = ['SPIRE500_err_EN1_croped.fits','SPIRE500_err_CDFS.fits']

list_names = ['all_1.38_mag_5']
map_names = ['SPIRE500_mask4K.fits']
emap_names = ['SPIRE500_err4.fits']

'''
list_names = ['r2']#,'all_1.77_mag']
map_names = ['SPIRE250_mask4.fits']
emap_names = ['SPIRE250_err4.fits']

#list_names = ['K_uds']#,'all_1.77_mag']
#map_names = ['SPIRE500_mask_uds.fits']
#emap_names = ['uds_cutout_noise_PLW.fits']


#list_names = ['H_huble'] #,K_uds']#,'all_1.77_mag']
#map_names = ['SPIRE500_mask_goods_croped.fits']
#emap_names = ['SPIRE500_err_goods_croped.fits']

##################################################
band = '250'
JK = 0
#PSFE = 19.0

print('JK',JK)

redshift_nodes  = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,10]) 
nuse_all = [0,1,2,3,4,5,6,7,8]

loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'
loc_cat = '/Users/Steven/Documents/python_simstack/cat/'


if band == '250':
  psf = 17.5 ## 17.5 #17.6
if band == '350':
  psf = 23.7 #23.9
if band == '500':
  psf = 34.6 #35.2

#psf = 5.0

fwhm = psf

#if PSFE > 17.5:
#  s1 = psf/(2*np.sqrt(2*np.log(2)))
#  s2 = PSFE/(2*np.sqrt(2*np.log(2)))

  #smL = np.sqrt(PSFE**2 - psf**2)
#  smL = np.sqrt(s2**2 - s1**2)
#  fwhm = PSFE

pixsize = 4.0 #3.0 #4.0
#kern = gauss_kern(psf, np.floor(psf * 18.)/pixsize, pixsize)
###len_ax = 30
len_ax = 15
#len_ax = 60

x, y = np.meshgrid(np.linspace(-1*len_ax*pixsize,len_ax*pixsize,2*len_ax+1), 
                   np.linspace(-1*len_ax*pixsize,len_ax*pixsize,2*len_ax+1))
d = np.sqrt(x*x+y*y)
print(fwhm)
sigma = fwhm/(2*np.sqrt(2*np.log(2)))
kern = np.exp(-( (d)**2 / ( 2.0 * sigma**2 ) ) )

#print(smL)
#sigma2 = smL #/(2*np.sqrt(2*np.log(2)))
#kernS = np.exp(-( (d)**2 / ( 2.0 * sigma2**2 ) ) )

zA = np.zeros(2)

#for nuse in nuse_all:
for q in range(0,np.size(list_names)):
#for q in range(3):



  if map_names[q]== emap_names[q]:
    emap, hds = fits.getdata(loc_m+emap_names[q],2, header = True)
    cmap, chd = fits.getdata(loc_m+map_names[q], 1, header = True)
  else:  
    emap, hds = fits.getdata(loc_m+emap_names[q],0, header = True)
    cmap, chd = fits.getdata(loc_m+map_names[q], 0, header = True)


  
  '''
  use = cmap > -999
  print(16.*np.size(cmap[use])/(3600**2))
  num_use = 1050
  num_use = 600
  #num_use = 550
  if JK == 1:
    cmap[:num_use,:num_use] = np.nan
  if JK == 2:
    cmap[num_use:,:num_use] = np.nan
  if JK == 3:
    cmap[:num_use,num_use:] = np.nan
  if JK == 4:
    cmap[num_use:,num_use:] = np.nan

  use = cmap > -999
  print(16.*np.size(cmap[use])/(3600**2))
  '''
  tbl = np.array(pd.read_table(loc_cat+list_names[q]+'.csv',sep=','))
  mag = tbl[:,3]
  ra_tot = tbl[:,1]
  dec_tot = tbl[:,2]
  z = tbl[:,5]
  
  '''
  if JK == 1:
    use_JK = (ra_tot > np.mean(ra_tot)) & (dec_tot > np.mean(dec_tot))
  if JK == 2:
    use_JK = (ra_tot > np.mean(ra_tot)) & (dec_tot <= np.mean(dec_tot))
  if JK == 3:
    use_JK = (ra_tot <= np.mean(ra_tot)) & (dec_tot > np.mean(dec_tot))
  if JK == 4:
    use_JK = (ra_tot <= np.mean(ra_tot)) & (dec_tot <= np.mean(dec_tot))
  mag = mag[~use_JK]
  ra_tot = ra_tot[~use_JK]
  dec_tot = dec_tot[~use_JK]
  z = z[~use_JK]
  '''

  radius = 1.1


  zeromask = np.ones(np.shape(cmap))
  ind_map_zero = np.where(clean_nans(cmap) == 0.0) ## in cmap all nan are 0
  zeromask[ind_map_zero] = 0.0

  
  #if PSFE > 17.5:
  #  print('in colcole', np.shape(cmap))
  #  cmap = scipy.signal.convolve2d(cmap, kernS, mode='same')
  #  print('in colcole', np.shape(cmap))

  cms = np.shape(cmap)    
  cw = WCS(chd, naxis = 2) 



  #aaa = mag < 23.0
  #ra_tot = ra_tot[aaa]
  #dec_tot = dec_tot[aaa]
  #z = z[aaa]
  #mag = mag[aaa]
  
  #print np.max(z), np.min(z)
  
  ################################## rm near 
  '''
  m_tot =mag
  dist_use = 4.0
  c_all = SkyCoord(ra=ra_tot*units.degree, dec=dec_tot*units.degree) 
  table = match_coordinates_sky(c_all, c_all, nthneighbor=2) 
  table2 = np.array(table)
  dist_arc = table2[1,:]*3600
  far = dist_arc > dist_use
  near = dist_arc <= dist_use

  ra_keep = ra_tot[far]
  dec_keep = dec_tot[far]
  m_keep = m_tot[far]

  ra_ap = ra_tot[near]
  dec_ap = dec_tot[near]
  mag_ap = m_tot[near]
  c_ap = c_all[near]

  more = True
  start = 0
  while (more == True) and (np.size(c_ap) != 0) :
    table = match_coordinates_sky(c_ap, c_ap, nthneighbor=2) 
    table2 = np.array(table)
    dist_arc = table2[1,:]*3600
    idx = table2[0,:].astype(int)
    good = (mag_ap < mag_ap[idx]) | (dist_arc > dist_use)

    ra_ap = ra_ap[good]
    dec_ap = dec_ap[good]
    mag_ap = mag_ap[good]
    c_ap = c_ap[good]
    if np.size(ra_ap) == start:
      more = False
    start = np.size(ra_ap)

  ra_tot = np.append(ra_keep, ra_ap)
  dec_tot = np.append(dec_keep, dec_ap)
  mag = np.append(m_keep, mag_ap)  
  '''
  ##################################
  
  #use = (m_tot > 0.000000001) & (z > redshift_nodes[nuse]) & (z <= redshift_nodes[nuse+1])
  #mag = -2.5*np.log10(m_tot[use]/3631)
  #ra_tot = ra_tot[use]
  #dec_tot = dec_tot[use]
  
  if list_names[q] == 'all_1.77_test':
    SOURCE_LIST = np.array([4,8,8.5,9,9.5,10,10.5,11,11.5,20])
  else:
    SOURCE_LIST = np.append(np.linspace(12.4,18,15),np.linspace(18.4,27.6,24))
    while np.max(mag) < SOURCE_LIST[np.size(SOURCE_LIST)-2]:
      SOURCE_LIST = SOURCE_LIST[0:-1]
    SOURCE_LIST[np.size(SOURCE_LIST)-1] = np.max(mag) + 1e-5
  '''  
  use = (z > redshift_nodes[nuse]) & (z <= redshift_nodes[nuse+1])
  ra_tot = ra_tot[use]
  dec_tot = dec_tot[use]
  mag = mag[use]
  print np.size(mag), np.size(ra_tot)
  '''

  run = np.size(SOURCE_LIST)
  nlists = run+2 - 1
  layers=np.zeros([nlists,cms[0],cms[1]])
  ngals_layer = np.zeros(nlists)


  for k in range(2,nlists):
      print k, 
      use = (mag >= SOURCE_LIST[k-2]) & (mag < SOURCE_LIST[k-2+1])
      ra = ra_tot[use]
      dec = dec_tot[use]
      if np.size(ra) > 0.5:
        ty,tx = cw.wcs_world2pix(ra, dec, 0)
        # CHECK FOR SOURCES THAT FALL OUTSIDE MAP
        ind_keep = np.where((np.round(tx) >= 0) & (np.round(tx) < cms[0]) & (np.round(ty) >= 0) & (np.round(ty) < cms[1]))
        real_x=np.round(tx[ind_keep]).astype(int)
        real_y=np.round(ty[ind_keep]).astype(int)
        # CHECK FOR SOURCES THAT FALL ON ZEROS
        ind_nz=np.where(cmap[real_x,real_y] != 0 )
      
        nt = np.shape(ind_nz)[1]
        ngals_layer[k] = nt
        print nt, np.size(tx)
        #print 'ngals: ' + str(nt)
        if nt > 0:
          real_x = real_x[ind_nz]
          real_y = real_y[ind_nz]
          for ni in range(nt):
            layers[k, real_x[ni],real_y[ni]]+=1.0

  
  DAT = np.array([np.append(zA,SOURCE_LIST[1:]),ngals_layer])
  ##np.savetxt(loc+'ns_'+list_names[q]+'.cat',DAT.T, delimiter="   ", fmt='%s', newline='\n')
  ##np.savetxt(loc+'ns_'+list_names[q]+'_JK'+str(JK)+'.cat',DAT.T, delimiter="   ", fmt='%s', newline='\n')
  
  ###np.savetxt(loc+'ns_'+list_names[q]+'_PSF'+str(PSFE)+'.cat',DAT.T, delimiter="   ", fmt='%s', newline='\n')

  #np.savetxt(loc+'ns_'+list_names[q]+'_FWHM'+str(psf)+'.cat',DAT.T, delimiter="   ", fmt='%s', newline='\n')
  #np.savetxt(loc+'ns_'+list_names[q]+'_'+str(redshift_nodes[nuse+1])+'.cat',DAT.T, delimiter="   ", fmt='%s', newline='\n')

  cfits_flat = np.asarray([])

  # STEP 2  - Convolve Layers and put in pixels`
  flattened_pixmap = np.sum(layers,axis=0)
  print('hier', fwhm,pixsize)
  total_circles_mask = circle_mask(flattened_pixmap, radius * fwhm, pixsize)
  ind_fit = np.where((total_circles_mask >= 1) & (zeromask != 0))
  del total_circles_mask
  nhits = np.shape(ind_fit)[1]
  cfits_flat = np.asarray([])
    
   
  once = True
  twice = True
  for u in range(nlists):
    print 'nsource', np.sum(layers[u,:,:])
    layer = layers[u,:,:]
    if u > 1.5:
      if np.sum(layer) > 0.5:
        tmap = scipy.signal.convolve2d(layer, kern, mode='same')
      else:
        tmap = layer  
    else:
      if (u == 0) :
        print u
        layer = layer + 1
        tmap = layer 
      elif (u == 1):
        if (np.size(ind_map_zero) > 0):
          print 'nanl'
          layer[ind_map_zero] = 1

        tmap = scipy.signal.convolve2d(layer, kern, mode='same')
    print np.min(tmap[ind_fit]) ,np.max(tmap[ind_fit])  

    print 'psf', np.sum(tmap[ind_fit]), np.sum(tmap), u
    cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))



  DAT = np.array([cmap[ind_fit],emap[ind_fit]]).T
  print np.shape(DAT)
  t = Table(DAT, names=('a','b'))
  #t.write(loc+'p_map_'+band+'_'+list_names[q]+'_MY'+str(redshift_nodes[nuse+1])+'.fits', format='fits')
  ##t.write(loc+'p_map_'+band+'_'+list_names[q]+'.fits', format='fits')
  #t.write(loc+'p_map_'+band+'_'+list_names[q]+'_JK'+str(JK)+'.fits', format='fits')

  #t.write(loc+'p_map_'+band+'_'+list_names[q]+'_PSF'+str(PSFE)+'.fits', format='fits')

  t.write(loc+'p_map_'+band+'_'+list_names[q]+'_FWHM'+str(psf)+'.fits', format='fits')

  t = Table([cfits_flat], names=('a'))
  #t.write(loc+'pm_'+band+'_beth_'+list_names[q]+'_MY'+str(redshift_nodes[nuse+1])+'.fits', format='fits')
  ##t.write(loc+'pm_'+band+'_beth_'+list_names[q]+'.fits', format='fits')
  t.write(loc+'pm_'+band+'_beth_'+list_names[q]+'_FWHM'+str(psf)+'.fits', format='fits')
  #t.write(loc+'pm_'+band+'_beth_'+list_names[q]+'_PSF'+str(PSFE)+'.fits', format='fits')
  #t.write(loc+'pm_'+band+'_beth_'+list_names[q]+'_JK'+str(JK)+'.fits', format='fits')




