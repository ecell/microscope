"""

tirfm_handler.py

"""
import sys
import os
import copy
import tempfile
import time
import math
import operator
import random
#import h5py

import scipy
import numpy

import parameter_configs
from epifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer
from effects_handler import PhysicalEffects

from scipy.special import j0
from scipy.misc    import toimage


class TIRFMConfigs(EPIFMConfigs) :

    '''

    TIRFM configration

	Evanescent Field
	    +
	Detector : EMCCD/CMOS
    '''

    def __init__(self, user_configs_dict = None):

        # default setting
        configs_dict = parameter_configs.__dict__.copy()
        #configs_dict_tirfm = tirfm_configs.__dict__.copy()
        #configs_dict.update(configs_dict_tirfm)

        # user setting
        if user_configs_dict is not None:
            if type(user_configs_dict) != type({}):
                print 'Illegal argument type for constructor of Configs class'
                sys.exit()
            configs_dict.update(user_configs_dict)

        for key, val in configs_dict.items():
            if key[0] != '_': # Data skip for private variables in setting_dict.
                if type(val) == type({}) or type(val) == type([]):
                    copy_val = copy.deepcopy(val)
                else:
                    copy_val = val
                setattr(self, key, copy_val)




    def set_Illumination_path(self) :

        #r = self.radial
        #d = self.depth
	r = numpy.linspace(0, 20000, 20001)
	d = numpy.linspace(0, 20000, 20001)

        # (plank const) * (speed of light) [joules meter]
        hc = 2.00e-25

	# Illumination : Assume that uniform illumination (No gaussian)
        w_0 = self.source_radius

	# power [joules/sec]
        P_0 = self.source_power

	# illumination area [m^2]
	A_0 = numpy.pi*w_0**2

        # single photon energy
        wave_length = self.source_wavelength*1e-9
        E_wl = hc/wave_length

	# photon flux [photons/sec]
        N_0 = P_0/E_wl

	################################################################
	# set for Evanescent field
        ################################################################

	# incident beam angle
	theta_in = self.source_angle

	sin2 = numpy.sin(theta_in)**2
	cos2 = numpy.cos(theta_in)**2

	# index refraction
	n_1 = 1.46 # fused silica
	n_2 = 1.33 # water (objective : water immersion)
	n2  = (n_2/n_1)**2

	# incident beam amplitude
	theta = numpy.pi/4.0
	A2_Is = N0/A_0*numpy.cos(theta)
	A2_Ip = N0/A_0*numpy.sin(theta)

	# Assume that the s-polarized direction is parallel to y-axis
	A2_z = A2_Ip*(4*cos2*(sin2 - n2)/(n2**2*cos2 + sin2 - n2))
	A2_y = A2_Is*(4*cos2/(1- n2))
	A2_x = A2_Ip*(4*cos2*sin2/(n2**2*cos2 + sin2 - n2))
	A2_T = A2_x + A2_y + A2_z

        I_r = numpy.array(map(lambda x : A2_T, r))

	# evanescent field depth (alpha = 1/depth)
	if (sin2/n > 1) : 
	    alpha = (4.0*numpy.pi/wave_length)*numpy.sqrt(sin2/n2 - 1)
	    I_d = numpy.exp(-alpha*d*1e-9)

	else :
	    I_d = numpy.array(map(lambda x : 1.00, d))

	# photon flux density [photon/(sec m^2)]
        self.source_flux = numpy.array(map(lambda x : I_r*x, I_d))

	print 'Photon Flux Density (Max) :', self.source_flux[0][0]



class TIRFMVisualizer(EPIFMVisualizer) :

	'''
	TIRFM visualization class
	'''

	def __init__(self, configs=TIRFMConfigs(), effects=PhysicalEffects()) :

		assert isinstance(configs, TIRFMConfigs)
		self.configs = configs

                assert isinstance(effects, PhysicalEffects)
                self.effects = effects

		"""
		Check and create the folder for image file.
		"""
		if not os.path.exists(self.configs.image_file_dir):
		    os.makedirs(self.configs.image_file_dir)
		#else:
		#    for file in os.listdir(self.configs.movie_image_file_dir):
		#	os.remove(os.path.join(self.configs.movie_image_file_dir, file))

                """
                set optical path from source to detector
                """
		self.configs.set_Optical_path()

