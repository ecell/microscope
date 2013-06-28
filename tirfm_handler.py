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
import h5py

#import pylab
import scipy
import numpy

import parameter_configs
#from epifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer
from pEpifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer
from effects_handler import PhysicalEffects

from scipy.special import j0
from scipy.misc    import toimage

#from matplotlib.backends.backend_pdf import PdfPages

IMAGE_SIZE_LIMIT=3000


class TIRFMConfigs(EPIFMConfigs) :

    '''
    TIRFM configration

	EPIFM configuration
	    +
	Mirror position (> critical angle )
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




    def set_EvanescentField(self, depth = None) :

        print '--- Evanescent Field :'
        print '\tPenetration Depth = ', depth, 'nm'

        self.penetration_depth = depth



    def set_Illumination_path(self) :

        r = self.radial
        z = self.depth

        # (plank const) * (speed of light) [joules meter]
        hc = 2.00e-25

        # (1) light source
        M2 = self.source_M2factor
        w_source = self.source_radius

        # power [joules/sec]
        P_0 = self.source_power

        # single photon energy
        wave_length = self.source_wavelength*1e-9
        E_wl = hc/wave_length

        # photon per sec [photons/sec]
        N_0 = P_0/E_wl

        # (2) beam expander
#        f_1 = self.expander_focal_length1
#        f_2 = self.expander_focal_length2

#        w_p = self.expander_pinhole_radius

#        w_BE = (f_2/f_1)*w_source

        # (3) scan and tube lens
#        f_s  = self.scanlens_focal_length
#        f_t1 = self.tubelens_focal_length1

#        w_tube = (f_t1/f_s)*w_BE

        # (4) objective
#        f_obj = self.objective_focal_length

        # Rayleigh range
#        z_R = numpy.pi*w_tube**2/wave_length

        # object distance to maximize image distance
#        s_obj = f_obj + z_R
#        w_obj = w_tube/numpy.sqrt((1 - s_obj/f_obj)**2 + (z_R/f_obj)**2)


        # (I) Beam Flux [photons/(m^2 sec)] (r <--> z)
#        w_r = w_obj*numpy.sqrt(1 + ((wave_length*r*1e-9)/(numpy.pi*w_obj**2))**2)
#        N_r = N_0*(1 - numpy.exp(-2*(w_p/w_r)**2))

        #bsf = numpy.array(map(lambda x, y : (2*x)/(numpy.pi*y**2)*numpy.exp(-2*(z*1e-9/y)**2), N_r, w_r))
        #flux = numpy.array(map(lambda x : (2*N_r)/(numpy.pi*w_r**2)*numpy.exp(-2*(x*1e-9/w_r)**2), z))
        flux = numpy.array(map(lambda x : (2*N_0)/(numpy.pi*w_source**2), z))

        # (II) Penetration Depth Function
#        n1 = self.objective_Ng
#        n2 = self.objective_Nm
#
#        length = math.sqrt(1 - self.objective_sin_max**2)/self.objective_sin_max
#        tan_th = self.mirror_position/length
#        sin2   = tan_th/math.sqrt(1 + tan_th**2)
#
#        ref = n1**2*sin2 - n2**2
#
#        if (ref > 0) :
#
#            penetration_depth = wave_length/(4*math.pi*math.sqrt(ref))/1e-9
#
#            print '--- TIRF Configuration : '
#
#        else :
#
#            penetration_depth = float('inf')
#
#            print '--- EPIF Configuration : '
#
#        print '\tPenetration Depth = ', penetration_depth, 'nm'

        func_r  = numpy.array(map(lambda x : 1.00, r))
        func_z  = numpy.exp(-z/self.penetration_depth)
        func_pd = numpy.array(map(lambda x : func_r*x, func_z))

        # Beam flux
        self.source_flux = flux*func_pd

	print 'Photon Flux (Surface) :', self.source_flux[0][0]



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
		if not os.path.exists(self.configs.movie_image_file_dir):
		    os.makedirs(self.configs.movie_image_file_dir)
		#else:
		#    for file in os.listdir(self.configs.movie_image_file_dir):
	#		os.remove(os.path.join(self.configs.movie_image_file_dir, file))

                """
                Image Size and Boundary
                """
                self.img_width  = int(self.configs.detector_image_size[0])
                self.img_height = int(self.configs.detector_image_size[1])

                if self.img_width > IMAGE_SIZE_LIMIT or self.img_height > IMAGE_SIZE_LIMIT :
                        raise VisualizerErrror('Image size is bigger than the limit size')

                """
                set optical path from source to detector
                """
		self.configs.set_Optical_path()

