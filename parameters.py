import pdb
import numpy as np
import ConfigParser
import os
import logging
import pprint
import astropy.cosmology as ac
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import Planck15, z_at_value
import astropy.units as u
#import NDpredict
from utils import string_is_true
#from astropy.cosmology import Planck15 as cosmo

def get_params(param_file_path):
    """
    Get parameter values and return them in the form of a dictionary.

    Parameters
    ----------
    param_file_path : str
        path to parameter file

    Returns
    -------
    params : dict
    """

    config = ConfigParser.SafeConfigParser()
    config.read(param_file_path)

    # Get "raw" dictionaries from `config` object
    #pdb.set_trace()
    raw_params = dict(config.items('general'))
    raw_cosmo_params = dict(config.items('cosmology'))
    raw_pop_params = dict(config.items('populations'))
    #try:
    #    raw_pop_params = dict(config.items('populations'))
    #except:
    #    pass
    try:
        raw_cut_params = dict(config.items('cuts'))
    except:
        pass
    raw_io_params = dict(config.items('io'))
    raw_binning_params = dict(config.items('binning'))
    raw_maps_to_stack_params = dict(config.items('maps_to_stack'))
    raw_map_path_params = dict(config.items('map_path'))
    raw_map_file_params = dict(config.items('map_file'))
    raw_noise_path_params = dict(config.items('map_path'))
    raw_noise_file_params = dict(config.items('noise_file'))
    raw_beams_params = dict(config.items('beams'))
    raw_color_correction_params = dict(config.items('color_correction'))
    raw_catalogs_params = dict(config.items('catalogs'))

    raw_io_params['param_file_path'] = os.path.abspath(param_file_path) # Store parameter file path

    # Convert "raw" config dictionary to "organized" dictionary `params`
    params = get_general_params(raw_params)
    params['io'] = get_io_parameters(raw_io_params)
    params['cosmo'] = get_cosmology_parameters(raw_cosmo_params)
    params['populations'] = get_population_parameters(raw_pop_params,params)
    #try:
    #    params['populations'] = get_population_parameters(raw_pop_params,params)
    #except:
    #    pass
    try:
        params['cuts'] = get_cut_parameters(raw_cut_params)
    except:
        pass
    params['map_files'] = get_maps_parameters(raw_maps_to_stack_params,raw_map_path_params,raw_map_file_params)
    params['noise_files'] = get_maps_parameters(raw_maps_to_stack_params,raw_noise_path_params,raw_noise_file_params)
    params['wavelength'] = get_wavelength_parameters(raw_maps_to_stack_params)
    params['psfs'] = get_beams_parameters(raw_maps_to_stack_params,raw_beams_params)
    params['color_correction'] = get_color_correction_parameters(raw_maps_to_stack_params,raw_color_correction_params)
    params['catalogs'] = get_catalogs_parameters(raw_catalogs_params)
    params['bins'] = get_binning_parameters(raw_binning_params)
    params['library_keys'] = params['map_files'].keys()

    params['z_err'] = raw_params['z_err']

    logging.info("---------- PARAMETER VALUES ----------")
    logging.info("======================================")
    logging.info("\n" + pprint.pformat(params, indent=4) + "\n")

    #pdb.set_trace()
    return params

def get_general_params(raw_params):
    params = {} # Initialize parameter dictionary

    # Catalog specific names for keys
    try:
        params['zkey'] = raw_params['zkey']
        params['mkey'] = raw_params['mkey']
    except:
        params['zkey'] = 'z_peak'
        params['mkey'] = 'lmass'
    try:
        params['ra_key'] = raw_params['ra_key']
        params['dec_key'] = raw_params['dec_key']
    except:
        params['ra_key'] = 'ra'
        params['dec_key'] = 'dec'
    try:
        params['uv_key'] = raw_params['uv_key']
        params['vj_key'] = raw_params['vj_key']
    except:
        params['uv_key'] = 'rf_U_V'
        params['vj_key'] = 'rf_V_J'


    # Have a floating background level instead of removing mean
    try:
        params['float_background'] = raw_params['float_background']
    except:
        params['float_background'] = False

    #pdb.set_trace()
    # Type of galaxy split.  Default is UVJ star-forming / quiescent
    try:
        params['galaxy_splitting_scheme'] = raw_params['classification_scheme']
    except KeyError:
        params['galaxy_splitting_scheme'] = 'sf-qt'
    try:
        params['save_bin_ids'] = string_is_true(raw_params['save_bin_ids'])
    except:
        params['save_bin_ids'] = True

    # If running bootstrap
    if string_is_true(raw_params['bootstrap'].split()[0]) == True:
        params['bootstrap'] = True
        params['boot0'] = float(raw_params['bootstrap'].split()[1])
        params['number_of_boots'] = float(raw_params['bootstrap'].split()[2])
        try:
            params['perturb_z'] = string_is_true(raw_params['bootstrap'].split()[3])
        except:
            params['perturb_z'] = False
        try:
            params['index_boots'] = string_is_true(raw_params['bootstrap'].split()[4])
        except:
            params['index_boots'] = False
        #pdb.set_trace()
    else:
        params['bootstrap'] = False
        params['boot0'] = 0
        params['number_of_boots'] = 1
        params['perturb_z'] = False

    return params

def get_wavelength_parameters(raw_maps_to_stack_params):
    wavelengths = {}

    for imap in raw_maps_to_stack_params:
        if string_is_true(raw_maps_to_stack_params[imap].split()[1]) == True:
            wavelengths[imap] = float(raw_maps_to_stack_params[imap].split()[0])

    return wavelengths

def get_binning_parameters(raw_params):
    binning = {}
    # Style of binning, optimal or evenly, and the number of bins (optional).
    # If number_of_bins not provided, will be decided by the binning code.
    #try:
    #    binning['optimal_binning'] = raw_params['optimal_binning'].split()[0]
    #except KeyError:
    #    binning['optimal_binning'] = False
    try:
        binning['optimal_binning'] = is_true(raw_params, 'optimal_binning')
    except KeyError:
        binning['optimal_binning'] = False

    #Optional number of bins
    if len(raw_params['optimal_binning'].split()) > 1:
        try:
            binning['number_of_bins'] = raw_params['optimal_binning'].split()[1]
        except KeyError:
            pass

    # If binning masses by Number Densities
    try:
        binning['bin_in_number_density'] = is_true(raw_params, 'bin_in_number_density')
    except KeyError:
        binning['bin_in_number_density'] = False

    # If binning redshifts in lookback time.
    try:
        binning['bin_in_lookback_time'] = is_true(raw_params, 'bin_in_lookback_time')
    except KeyError:
        binning['bin_in_lookback_time'] = False

    # If stacking entire catalog at once, rather that in redshift slices.
    # Still unclear if this is advantageous or not.
    try:
        binning['stack_all_z_at_once'] = is_true(raw_params, 'all_z_at_once')
    except KeyError:
        binning['stack_all_z_at_once'] = False

    if is_true(raw_params,'optimal_binning') == False:
        z_nodes = []
        m_nodes = []
        for i in raw_params['redshift_nodes'].split():
            z_nodes.append(float(i))
        for j in raw_params['mass_nodes'].split():
            m_nodes.append(float(j))
        if binning['bin_in_number_density'] == True:
            nd_nodes = []
            try:
                for j in raw_params['number_density_nodes'].split():
                    nd_nodes.append(float(j))
            except:
                for j in raw_params['mass_nodes'].split():
                    nd_nodes.append(float(j))


        binning['t_nodes'] = z_nodes
        binning['z_nodes'] = z_nodes
        binning['m_nodes'] = m_nodes
        #binning['nd_nodes'] = nd_nodes

        # This is tough...
        if binning['bin_in_number_density'] == True:
            pdb.set_trace()
            binning['m_nodes'] = np.array([NDpredict.getmass_illustris(nd_nodes,z_mid[i]) for i in m_nodes])

        if binning['bin_in_lookback_time'] == True:
            binning['t_nodes'] = z_nodes
            binning['z_nodes'] = np.array([z_at_value(Planck15.age,(cosmo.age(0).value - i) * u.Gyr) for i in z_nodes])

    return binning

def get_io_parameters(raw_params):
    io = {}

    io['param_file_path']            = raw_params['param_file_path']
    try:
        io['shortname']              = raw_params['shortname']
    except KeyError:
        io['shortname']              = ''

    io['output_folder']              = os.environ[raw_params['output_folder'].split()[0]] + raw_params['output_folder'].split()[1]
    io['flux_densities_filename']    = raw_params['flux_densities_filename']
    #io['flux_densities_filename']    = 'simstack_flux_densities'

    return io

def get_cosmology_parameters(raw_params):
    '''
    Returns
    -------
    cosmo : astropy.cosmology object
        object containing cosmological parameters
    '''
    omega_m0    = float(raw_params['omega_m'])    # Present-day matter density
    omega_l0    = float(raw_params['omega_l'])    # Present-day dark energy density
    omega_k0    = float(raw_params['omega_k'])    # Present-day spatial curvature density
    hubble_h0   = float(raw_params['h'])          # Present-day reduced Hubble constant: h0 = H0 / (100 km/s/Mpc)

    H0          = hubble_h0*100.
    cosmo       = ac.LambdaCDM(H0=H0, Om0=omega_m0, Ode0=omega_l0)

    return cosmo

def get_cut_parameters(raw_cut_params):
    cuts_dict = {}
    # Special case for 5pop and 4pop
    for pop in raw_cut_params:
        cuts_dict[pop] = float(raw_cut_params[pop])
    return cuts_dict

def get_population_parameters(raw_pop_params, params):

    cuts_dict = {}

    if params['galaxy_splitting_scheme'] == 'general':
        for pop in raw_pop_params:
            tst = [int(raw_pop_params[pop][0])]
            if len(raw_pop_params[pop].split()) > 1:
                tst.append([k for k in raw_pop_params[pop][1:].split()])
                for k in range(len(tst[1])):
                    try:
                        bl = string_is_true(tst[1][k])
                        tst[1][k] = bl
                        #print 'is a boolean'
                    except NameError:
                        try:
                            float(tst[1][1])
                            tst[1][k]=float(tst[1][k])
                            #print 'is a float'
                        except ValueError:
                            #print 'do nothin'
                            pass
            else:
                tst.append([])

            cuts_dict[pop] = tst

    elif params['galaxy_splitting_scheme'] == 'sf-qt' or params['galaxy_splitting_scheme'] == '4pops' or params['galaxy_splitting_scheme'] == '5pops':
        for pop in raw_pop_params:
            cuts_dict[pop] = float(raw_pop_params[pop])

    elif params['galaxy_splitting_scheme'] == 'uvj':
        cuts_dict['c_nodes'] = [float(n) for n in raw_pop_params['uvj_nodes'].split()]
        cuts_dict['pop_names'] = [n for n in raw_pop_params['pop_names'].split()]

    return cuts_dict

def get_maps_parameters(raw_maps_to_stack_params,raw_map_path_params,raw_map_file_params):
    maps = {}

    for imap in raw_maps_to_stack_params:
        if string_is_true(raw_maps_to_stack_params[imap].split()[1]) == True:
            maps[imap] = os.environ[raw_map_path_params[imap].split()[0]] + raw_map_path_params[imap].split()[1] + raw_map_file_params[imap]

    return maps

def get_beams_parameters(raw_maps_to_stack_params,raw_beams_params):
    psfs = {}

    for imap in raw_maps_to_stack_params:
        if string_is_true(raw_maps_to_stack_params[imap].split()[1]) == True:
            psfs[imap+'_beam_area'] = float(raw_beams_params[imap].split()[1])
            if is_float(raw_beams_params[imap].split()[0]) == True:
                psfs[imap+'_fwhm'] = float(raw_beams_params[imap].split()[0])
            else:
                psfs[imap+'_beam_file'] = raw_beams_params[imap].split()[0]

    return psfs

def get_color_correction_parameters(raw_maps_to_stack_params,raw_color_correction_params):
    color_correction = {}

    for imap in raw_maps_to_stack_params:
        if string_is_true(raw_maps_to_stack_params[imap].split()[1]) == True:
            color_correction[imap+''] = float(raw_color_correction_params[imap])

    return color_correction

def get_catalogs_parameters(raw_catalog_params):
    catalog = {}
    try:
        catalog['catalog_path'] = os.environ[raw_catalog_params['catalog_path'].split()[0]] + raw_catalog_params['catalog_path'].split()[1]
    except:
        catalog['catalog_path'] = raw_catalog_params['catalog_path']

    catalog['catalog_file'] = raw_catalog_params['catalog_file']
    if 'features_file' in raw_catalog_params:
        catalog['features_file'] = raw_catalog_params['features_file']

    return catalog

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def string_is_true(sraw):
    """Is string true? Returns boolean value.
    """
    s       = sraw.lower() # Make case-insensitive

    # Lists of acceptable 'True' and 'False' strings
    true_strings    = ['true', 't', 'yes', 'y', '1']
    false_strings    = ['false', 'f', 'no', 'n', '0']
    if s in true_strings:
        return True
    elif s in false_strings:
        return False
    else:
        logging.warning("Input not recognized for parameter: %s" % (key))
        logging.warning("You provided: %s" % (sraw))
        raise

def is_true(raw_params, key):
    """Is raw_params[key] true? Returns boolean value.
    """
    sraw    = raw_params[key]
    s       = sraw.lower() # Make case-insensitive

    # Lists of acceptable 'True' and 'False' strings
    true_strings    = ['true', 't', 'yes', 'y', '1']
    false_strings    = ['false', 'f', 'no', 'n', '0']
    if s in true_strings:
        return True
    elif s in false_strings:
        return False
    else:
        logging.warning("Input not recognized for parameter: %s" % (key))
        logging.warning("You provided: %s" % (sraw))
        raise

### FOR TESTING ###
if __name__=='__main__':
    import os, sys
    import pprint

    param_fp = sys.argv[1]
    print("")
    print("Testing %s on %s..." % (os.path.basename(__file__), param_fp))
    print("")
    pprint.pprint(get_params(param_fp))
