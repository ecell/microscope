"""

fcs_handler.py

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


class FCSConfigs(EPIFMConfigs) :

    '''
    FCS configuration

	EPIFM configuration
	    +
	Pinhole Lens
	Detector : PMT, ADP ...etc
    '''

    def __init__(self, user_configs_dict = None):

        # default setting
        configs_dict = parameter_configs.__dict__.copy()
        #configs_dict_fcs = fcs_configs.__dict__.copy()
        #configs_dict.update(configs_dict_fcs)

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



    def set_PinholeLens(self, radius = None,
                        focal_length = None) :

        print '--- Pinhole Lens :'

        self._set_data('pinholelens_switch', True)
        self._set_data('pinholelens_radius', radius)
        self._set_data('pinholelens_focal_length', focal_length)

        print '\tPinhole Radius = ', self.pinholelens_radius, 'm'
        print '\tFocal Length = ', self.pinholelens_focal_length, 'm'




    def set_Detection_path(self) :

        wave_length = self.psf_wavelength*1e-9

        # source waist aat image position
        waist = self.source_divergence[0]

        # Magnification
        Mag = 1.0

        # (2) objective lens
        f_obj = self.objective_focal_length
        w_obj = waist*numpy.sqrt(1 + ((wave_length*f_obj)/(numpy.pi*waist**2))**2)

        # (3) tube lens
        f_t = self.tubelens_focal_length
        w_t = (wave_length*f_obj)/(numpy.pi*w_obj)

        # Magnification : Obj to tube lens
        Mag = (f_t/f_obj)*Mag

        # (4) scan lens
        f_s = self.scanlens_focal_length
        w_s = w_t*numpy.sqrt(1 + ((wave_length*f_t)/(numpy.pi*w_t**2))**2)

        # (5) pinhole lens in front of detector
        f_p = self.pinhole_focal_length
        w_p = (wave_length*f_s)/(numpy.pi*w_s)

        # Magnification : scan to pinhole lens
        Mag = (f_p/f_s)*Mag

        # set image scaling factor
        voxel_radius = self.spatiocyte_VoxelRadius

        view = (2.0*self.pinholelens_radius)/(2.0*voxel_radius)
        zoom = self.detector_zoom

        self.image_scaling = view/(Mag*zoom)

        # Detector PSF
        self.set_PSF_detector(Mag)




class FCSVisualizer(EPIFMVisualizer) :

	'''
	FCS Visualization class of e-cell simulator
	'''

	def __init__(self, configs=FCSConfigs(), output_file=None) :

		assert isinstance(configs, FCSConfigs)
		self.configs = configs
		self.output_file = output_file

		"""
		Check and create the folder for image file.
		"""
		if not os.path.exists(self.configs.movie_image_file_dir):
		    os.makedirs(self.configs.movie_image_file_dir)
		else:
		    for file in os.listdir(self.configs.movie_image_file_dir):
			os.remove(os.path.join(self.configs.movie_image_file_dir, file))

                """
		Output Data
                """
                with open(self.output_file, 'w') as output :
                    output.write('#time\tintensity\t\n')
                    output.write('\n')

                """
                Optical Path
                """
		self.configs.set_Optical_path()



        def output_frames(self, num_div=1) :

                # define observational image plane in nm-scale
                voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9

                Nz = int(self.configs.spatiocyte_lengths[0][2]*voxel_size)
                Ny = int(self.configs.spatiocyte_lengths[0][1]*voxel_size)
                Nx = int(self.configs.spatiocyte_lengths[0][0]*voxel_size)

                z = numpy.linspace(0, Nz-1, Nz)
                y = numpy.linspace(0, Ny-1, Ny)
                Z, Y = numpy.meshgrid(z, y)

                # focal point
                p_0 = numpy.array([Nx, Ny, Nz])*self.configs.detector_focal_point

                # beam position : Assuming beam position = focal point for temporary
                p_b = copy.copy(p_0)

		# frame interval
		frame_interval = 1.0/self.configs.detector_fps
		exposure_time  = self.configs.detector_exposure_time

		time = self.configs.detector_start_time
		end  = self.configs.detector_end_time

		# data-time interval
                t0 = self.configs.spatiocyte_data[0][0]
                t1 = self.configs.spatiocyte_data[1][0]

		delta_data  = t1 - t0
                delta_frame = int(frame_interval/delta_data)

                # create frame data composed by frame element data
		count = 0

		while (time < end) :

                    # set image file name
                    self.image_file_name = os.path.join(self.configs.movie_image_file_dir, \
                                                self.configs.movie_image_file_name_format % (count))

                    print 'time : ', time, ' sec'

                    cell = numpy.zeros(shape=(Nz, Ny))

		    count_start = count*delta_frame
		    count_end = (count + 1)*delta_frame

		    frame_data = self.configs.spatiocyte_data[count_start:count_end]

		    # loop for frame data
                    for i in range(len(frame_data)) :

			# i-th data in a frame
			i_time = frame_data[i][0]
			data   = frame_data[i][1]
			total  = len(data)

			Norm = numpy.exp(-(i_time - time)/exposure_time)

			# loop for particles
			for j in range(total) :

			    c_id, s_id, l_id = data[j]

                            # particles coordinate in real(nm) scale
                            pos = self.get_coordinate(c_id)
                            p_i = numpy.array(pos)*voxel_size

                            print j, p_i
                            # get signal matrix
                            signal = Norm*numpy.array(self.get_signal(p_i, p_b, p_0))

                            # add signal matrix to image plane
                            cell = self.overwrite(cell, signal, p_i)

		    # intensity
                    intensity = self.detector_output(cell)

                    # Detector Output
                    data_line  = str(time) + '\t'
                    data_line += str(intensity) + '\t'
                    data_line += '\n'

                    with open(self.output_file, 'a') as output :
                        output.write(data_line)


		    time  += frame_interval
		    count += 1



	def detector_output(self, cell) :

		# Detector Output
                voxel_size = (2.0*self.configs.spatiocyte_VoxelRadius)/1e-9

		Nw_pixel = 1#self.img_width
		Nh_pixel = 1#self.img_height

		Np = int(self.configs.image_scaling*voxel_size)

                # image in nm-scale
		Nw_camera = Nw_pixel*Np
		Nh_camera = Nh_pixel*Np

		Nw_cell = len(cell)
		Nh_cell = len(cell[0])

		if (Nw_camera > Nw_cell) :

			w_cam_from = int((Nw_camera - Nw_cell)/2.0)
			w_cam_to   = w_cam_from + Nw_cell
                        w_cel_from = 0
                        w_cel_to   = Nw_cell

		else :

                        w_cam_from = 0
                        w_cam_to   = Nw_camera
                        w_cel_from = int((Nw_cell - Nw_camera)/2.0)
                        w_cel_to   = w_cel_from + Nw_camera

		if (Nh_camera > Nh_cell) :

                        h_cam_from = int((Nh_camera - Nh_cell)/2.0)
                        h_cam_to   = h_cam_from + Nh_cell
                        h_cel_from = 0
                        h_cel_to   = Nh_cell

                else :

                        h_cam_from = 0
                        h_cam_to   = int(Nh_camera)
                        h_cel_from = int((Nh_cell - Nh_camera)/2.0)
                        h_cel_to   = h_cel_from + Nh_camera


		# image in nm-scale
		plane = cell[w_cel_from:w_cel_to, h_cel_from:h_cel_to]

		# convert image in nm-scale to pixel-scale
		#signal = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])
		signal = numpy.sum(plane)
		intensity = self.A2D_converter(signal)

		return intensity



