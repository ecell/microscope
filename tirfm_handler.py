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

import pylab
import scipy
import numpy

import parameter_configs
from epifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer

from scipy.special import j0
from scipy.misc    import toimage

from matplotlib.backends.backend_pdf import PdfPages

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




    def set_Mirror(self, position = None ) :

        print '--- Mirror :'

        self._set_data('mirror_position', position)

        print '\tMirror Position = ', self.mirror_position


    def set_Illumination_path(self) :

        r = self.radial
        z = self.depth

        # (plank const) * (speed of light) [joules meter]
        hc = 2.00e-25

        # (1) light source
        M2  = self.source_M2factor
        width_source = self.source_radius
        P_0 = self.source_power
        wave_length = self.source_wavelength*1e-9

        # calculate for single photon energy
        E_wl = hc/wave_length

        # convertion factor to the photon number
        voxel_area = numpy.pi*(self.spatiocyte_VoxelRadius**2)
        F0 = (voxel_area/E_wl)*self.detector_exposure_time

        # (2) beam expander
        f_1 = self.expander_focal_length1
        f_2 = self.expander_focal_length2

        w_1 = M2*(wave_length*f_1)/(numpy.pi*width_source)
        w_p = self.expander_pinhole_radius
        w_2 = w_1*numpy.sqrt(1 + ((wave_length*f_2)/(numpy.pi*w_1**2))**2)

        # (3) scan lens
        f_s = self.scanlens_focal_length
        w_s = (wave_length*f_s)/(numpy.pi*w_2)

        # (4) tube lens
        f_t = self.tubelens_focal_length
        w_t = w_s*numpy.sqrt(1 + ((wave_length*f_t)/(numpy.pi*w_s**2))**2)

        # (5) objective
        f_obj = self.objective_focal_length
        w_obj = (wave_length*f_obj)/(numpy.pi*w_t)


        # Illumination PSF : (I) Beam Spread Function
        w_z = w_obj*numpy.sqrt(1 + ((wave_length*z*1e-9)/(numpy.pi*w_obj**2))**2)
        P_z = P_0*(1 - numpy.exp(-2*(w_p/w_z)**2))

        bsf = numpy.array(map(lambda x, y : (2*x)/(numpy.pi*y**2)*numpy.exp(-2*(r*1e-9/y)**2), P_z, w_z))

        self.source_divergence = w_z

        # Illumination PSF : (II) Penetration Depth Function
        n1 = self.objective_Ng
        n2 = self.objective_Nm

        length = math.sqrt(1 - self.objective_sin_max**2)/self.objective_sin_max
        tan_th = self.mirror_position/length
        sin2   = tan_th/math.sqrt(1 + tan_th**2)

        ref = n1**2*sin2 - n2**2

        if (ref > 0) :

            self.penetration_depth = wave_length/(4*math.pi*math.sqrt(ref))/1e-9

            print '--- TIRF Configuration : '

        else :

            self.penetration_depth = float('inf')

            print '--- EPIF Configuration : '

        print '\tPenetration Depth = ', self.penetration_depth, 'nm'

        psf_r  = numpy.array(map(lambda x : 1.00, r))
        psf_z  = numpy.exp(-z/self.penetration_depth)
        psf_pd = numpy.array(map(lambda x : psf_r*x, psf_z))

        # Illumination PSF : Total
        self.source_psf = F0*bsf*psf_pd




class TIRFMVisualizer(EPIFMVisualizer) :

	'''
	TIRFM visualization class
	'''

	def __init__(self, configs=TIRFMConfigs()) :

		assert isinstance(configs, TIRFMConfigs)
		self.configs = configs

		"""
		Check and create the folder for image file.
		"""
		if not os.path.exists(self.configs.movie_image_file_dir):
		    os.makedirs(self.configs.movie_image_file_dir)
		else:
		    for file in os.listdir(self.configs.movie_image_file_dir):
			os.remove(os.path.join(self.configs.movie_image_file_dir, file))

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



