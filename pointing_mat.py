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
from astropy import units as u
from astropy.coordinates import match_coordinates_sky

from sklearn.utils import shuffle



band = '250'

name2 = 'full_fwhm1'
name3 = 'full_noise_2arc_30arc'
name2 = 'full_noise_2arc_30arc'

#name2 = 'full_noise'

loc_m = '/Users/Steven/Documents/python_simstack/maps/'
map_250 = 'SPIRE'+band+'_gaussbeam_MINE_'+name2+'.fits'
map_ex = 'SPIRE'+band+'_gaussbeam_MINE_'+name3+'.fits'


if band == '250':
  psf = 17.6
if band == '350':
  psf = 23.9
if band == '500':
  psf = 35.2


cmap, chd = fits.getdata(loc_m+map_250, 0, header = True)
cmap_ex, chd_ex = fits.getdata(loc_m+map_ex, 0, header = True)
clean_nans(cmap_ex) == 0.0 ## in cmap all nan are 0


cms = np.shape(cmap)    
cw = WCS(chd, naxis = 2)

emap, hds = fits.getdata(loc_m+'SPIRE'+band+'_gaussbeam_MINE_noise.fits',0, header = True)
#emap, hds = fits.getdata(loc_m+'SPIRE250_gaussbeam_MINE_noise_tiny01.fits',0, header = True)

loc = '/Users/Steven/Documents/prepare_simsack/point_mat/'


fwhm = psf
pixsize = 4.0
#kern = gauss_kern(psf, np.floor(psf * 18.)/pixsize, pixsize)
len_ax = 30
x, y = np.meshgrid(np.linspace(-1*len_ax*pixsize,len_ax*pixsize,2*len_ax+1), 
                   np.linspace(-1*len_ax*pixsize,len_ax*pixsize,2*len_ax+1))
d = np.sqrt(x*x+y*y)
sigma = fwhm/(2*np.sqrt(2*np.log(2)))
kern = np.exp(-( (d)**2 / ( 2.0 * sigma**2 ) ) )

radius = 1.1

zeromask = np.ones(np.shape(cmap))
ind_map_zero = np.where(clean_nans(cmap) == 0.0) ## in cmap all nan are 0
zeromask[ind_map_zero] = 0.0



name = 'Mock_cat_Bethermin2017.fits'
loc3 = '/Users/Steven/Downloads/'
hdu = fits.open(loc3+name)
f2 = hdu[1].data
flux2 = (f2[0])[14]
#print np.sum(flux2)
m_tot = (f2[0])[10]
ra_tot = (f2[0])[1]
dec_tot = (f2[0])[2]

redshift_nodes  = np.array([0.00, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,10]) 
nuse =9

z =  (f2[0])[0]

###ra_tot, dec_tot = shuffle(ra_tot,dec_tot)


mat = np.loadtxt('/Users/Steven/Documents/prepare_simsack/pois_radec.cat')
#ra_tot =  mat[:,0]
#dec_tot = mat[:,1]


################################## rm near 
dist_use = 2.0
c_all = SkyCoord(ra=ra_tot*u.degree, dec=dec_tot*u.degree) 
table = match_coordinates_sky(c_all, c_all, nthneighbor=2) 
table2 = np.array(table)
dist_arc = table2[1,:]*3600
far = dist_arc > dist_use
near = dist_arc <= dist_use

ra_keep = ra_tot[far]
dec_keep = dec_tot[far]
m_keep = m_tot[far]
z_keep = z[far]

ra_ap = ra_tot[near]
dec_ap = dec_tot[near]
mag_ap = m_tot[near]
z_ap = z[near]
c_ap = c_all[near]

more = True
start = 0
while more == True:
  table = match_coordinates_sky(c_ap, c_ap, nthneighbor=2) 
  table2 = np.array(table)
  dist_arc = table2[1,:]*3600
  idx = table2[0,:].astype(int)
  good = (mag_ap > mag_ap[idx]) | (dist_arc > dist_use)

  ra_ap = ra_ap[good]
  dec_ap = dec_ap[good]
  mag_ap = mag_ap[good]
  z_ap = z_ap[good]
  c_ap = c_ap[good]
  if np.size(ra_ap) == start:
    more = False
  start = np.size(ra_ap)

ra_tot = np.append(ra_keep, ra_ap)
dec_tot = np.append(dec_keep, dec_ap)
m_tot = np.append(m_keep, mag_ap)
z = np.append(z_keep,z_ap)

##################################




#use = m_tot > 0.0000001
use = (m_tot > 0.000000001) & (z > redshift_nodes[nuse]) & (z <= redshift_nodes[nuse+1])

mag = -2.5*np.log10(m_tot[use]/3631)
#print np.sum(flux2[use])
#exit()
ra_tot = ra_tot[use]
dec_tot = dec_tot[use]




SOURCE_LIST = np.append(np.linspace(12.4,18,15),np.linspace(18.4,26.4,21))


#SOURCE_LIST = np.append(np.linspace(12.4,18,15),np.linspace(18.4,27.6,24))

#SOURCE_LIST = np.linspace(12.5,26.5,29)
#SOURCE_LIST = np.linspace(12.5,26.5,15)
#SOURCE_LIST = np.linspace(12.4,26.4,71)
#SOURCE_LIST = np.linspace(12.5,26.5,8)
#SOURCE_LIST = np.linspace(12.4,26.4,141)
#SOURCE_LIST = np.linspace(12.4,26.4,281)
#SOURCE_LIST = np.linspace(12.4,26.4,1401)
run = np.size(SOURCE_LIST)
nlists = run+2 - 1
layers=np.zeros([nlists,cms[0],cms[1]])
ngals_layer = np.zeros(nlists)




for k in range(2,nlists):
#for k in range(nlists-13,nlists):
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
      ###ind_nz=np.where(cmap[real_x,real_y] != 0 )
      ind_nz=np.where(cmap_ex[real_x,real_y] != 0 )

      #ind_nz=np.where(cmap[real_x,real_y] != -99 )
    
      nt = np.shape(ind_nz)[1]
      ngals_layer[k] = nt
      print nt, np.size(tx)
      #print 'ngals: ' + str(nt)
      if nt > 0:
        real_x = real_x[ind_nz]
        real_y = real_y[ind_nz]
        for ni in range(nt):
    	  layers[k, real_x[ni],real_y[ni]]+=1.0


#np.savetxt(loc+'ns_beth_2arc_30arc_2ac.cat',ngals_layer.T, delimiter="   ", fmt='%s', newline='\n')
np.savetxt(loc+'ns_beth_2arc_30arc_2ac_'+str(redshift_nodes[nuse+1])+'_2ac.cat',ngals_layer.T, delimiter="   ", fmt='%s', newline='\n')

cfits_flat = np.asarray([])

# STEP 2  - Convolve Layers and put in pixels`
flattened_pixmap = np.sum(layers,axis=0)
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
    print u
    ###tmap = smooth_psf(layer, kern)
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
      #tmap = smooth_psf(layer, kern)
  print np.min(tmap[ind_fit]) ,np.max(tmap[ind_fit])  

  print 'psf', np.sum(tmap[ind_fit]), np.sum(tmap), u
  cfits_flat = np.append(cfits_flat,np.ndarray.flatten(tmap[ind_fit]))
  #if  (u > 31) and (u < 33):
    #plt.hist(tmap[ind_fit])
    #plt.yscale('log')
    #aaa = tmap < 0.00001
    #tmap[aaa] =  1000
    #plt.imshow(tmap)
    #plt.show()


DAT = np.array([cmap[ind_fit],emap[ind_fit]]).T
print np.shape(DAT)
t = Table(DAT, names=('a','b'))
t.write(loc+'p_map_'+band+'_beth_'+name2+'2ac_MY'+str(redshift_nodes[nuse+1])+'.fits', format='fits')

t = Table([cfits_flat], names=('a'))
t.write(loc+'pm_'+band+'_beth_'+name2+'2ac_MY'+str(redshift_nodes[nuse+1])+'.fits', format='fits')





