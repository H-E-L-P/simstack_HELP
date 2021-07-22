# Standard modules
import pdb
import os
import os.path
import sys
import six
import shutil
import time
import logging
import importlib
import cPickle as pickle
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.wcs import WCS
import parameters
from utils import clean_args

class Bootstrap:

	def __init__(self,tbl):

		self.table = tbl

	def perturb_catalog(self, z_err, zkey, perturb_z = True, draw_z = False):
		'''
		Default is a simple bootstrap with replacement.
		If perturb_z is True, then the redshifts (z_peak) are perturbed by random normal with sigma z_err
		If draw_z is True, then the redshifts are drawn from the redshift PDF (output by EAZY) of each object (TODO!)
		'''
		pseudo_cat = self.table.copy()
		ngals = len(pseudo_cat)
		perturb_z = True
		if perturb_z == True:
			print 'start perturbed z'
			pseudo_z = self.table[zkey] + self.table[z_err]*np.random.randn(len(self.table[z_err]))
			#pseudo_z = self.table['z_peak'] + self.table['z_peak'].iloc[np.random.randint(0, ngals, size=ngals)]
			#self.table['z_err']*np.random.randn(len(self.table['z_err']))
			pseudo_cat[zkey] = pseudo_z
		#print pseudo_cat.sample(ngals,replace=True)
		#Simple Bootstrap.  Sample draws from pseudo_cat with replacement.
		self.pseudo_cat = pseudo_cat.sample(ngals,replace=False) # is true
