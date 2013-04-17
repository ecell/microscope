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
import csv

import pylab
import scipy
import numpy

import parameter_configs
#from epifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer
from pEpifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer

from scipy.special import j0
from scipy.misc    import toimage

from matplotlib.backends.backend_pdf import PdfPages

IMAGE_SIZE_LIMIT=3000


class KinesinConfigs(EPIFMConfigs) :

    '''
    Kinesin configration

	EPIFM configuration
	    +
	read off-lattice coordinate
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



    def set_DataFile(self, csv_file_path_list) :

        # read hdf5 lattice file
        for csv_file_path in csv_file_path_list :

            try :

                csv_file = open(csv_file_path, 'r')

		dataset = []

		for row in csv.reader(csv_file) :
		    dataset.append(row)

                ### header
		header = dataset[0]
		interval = float(header[0].split('=')[1])
		lengths  = [(float(header[3].split('=')[1]), float(header[2].split('=')[1]), float(header[1].split('=')[1]))]
		voxel_r  = float(header[4].split('=')[1])
		s_id = 6
		l_id = 0

                ### particle data in time-series
                data = []

		time = 0
                for i in range(1, len(dataset)) :

                    data_i = dataset[i]
		    particles = []

		    for j in range(1, len(data_i), 3) :

			c_id = (float(data_i[j]), float(data_i[j+1]), float(data_i[j+2]))
			particles.append((c_id, s_id, l_id))
			

		    element = [time, particles]
                    data.append(element)

		    time += interval


                data.sort(lambda x, y:cmp(x[0], y[0]))

                # get data
                self._set_data('spatiocyte_data', data)

                # get species properties
                #self._set_data('spatiocyte_species_id', map(lambda x : x[0], species))
                #self._set_data('spatiocyte_index',      map(lambda x : x[1], species))
                #self._set_data('spatiocyte_diffusion',  map(lambda x : x[3], species))
                #self._set_data('spatiocyte_radius',     map(lambda x : x[2], species))

                # get lattice properties
                #self._set_data('spatiocyte_lattice_id', map(lambda x : x[0], lattice))
                self._set_data('spatiocyte_lengths', lengths)
                self._set_data('spatiocyte_VoxelRadius', voxel_r)
                self._set_data('spatiocyte_theNormalizedVoxelRadius', 0.5)
                #self._set_data('spatiocyte_theStartCoord', lattice[0][4])
                #self._set_data('spatiocyte_theRowSize',    lattice[0][6])
                #self._set_data('spatiocyte_theLayerSize',  lattice[0][5])
                #self._set_data('spatiocyte_theColSize',    lattice[0][7])


            except Exception, e :
                        if not self.ignore_open_errors:
                            raise
                        print 'Ignoring error: ', e






class KinesinVisualizer(EPIFMVisualizer) :

	'''
	Kinesin visualization class
	'''

	def __init__(self, configs=KinesinConfigs()) :

		assert isinstance(configs, KinesinConfigs)
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


        def get_coordinate(self, aCoord) :

                """
                get (x, y, z) coordinate
                """
                point_y = aCoord[1]
                point_z = aCoord[2]
                point_x = aCoord[0]

                return point_x, point_y, point_z
                #return point_y, point_x, point_z

