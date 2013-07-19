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
import ctypes
import multiprocessing

#import pylab
import scipy
import numpy

import parameter_configs
#from epifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer
from pEpifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer

from time import sleep
from scipy.special import j0
from scipy.misc    import toimage

#from matplotlib.backends.backend_pdf import PdfPages

IMAGE_SIZE_LIMIT=3000


class ConfocalConfigs(EPIFMConfigs) :

    '''
    Confocal configuration

	EPIFM configuration
	    +
	Pinhole Lens
	Detector : PMT, ADP ... etc
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


    def set_Illumination_path(self) :

        r = self.radial
        z = self.depth

        # (plank const) * (speed of light) [joules meter]
        hc = 2.00e-25

        # (1) light source
        M2  = self.source_M2factor
        w_source = self.source_radius

        # power [joules/sec]
        P_0 = self.source_power

        # single photon energy
        wave_length = self.source_wavelength*1e-9
        E_wl = hc/wave_length

        # photon per sec [photons/sec]
        N_0 = P_0/E_wl

        # (2) beam expander
        f_1 = self.expander_focal_length1
        f_2 = self.expander_focal_length2

        w_p = self.expander_pinhole_radius

        w_BE = (f_2/f_1)*w_source

        # (3) scan and tube lens
        f_s = self.scanlens_focal_length
        f_t = self.tubelens_focal_length1

        w_tube = (f_t/f_s)*w_BE

        # (4) objective
        f_obj = self.objective_focal_length

        # Rayleigh range
        z_R = numpy.pi*w_tube**2/wave_length

        # object distance to maximize image distance
        s_obj = f_obj + z_R
        w_obj = w_tube/numpy.sqrt((1 - s_obj/f_obj)**2 + (z_R/f_obj)**2)

        # Beam Flux [photons/(m^2 sec)]
        w_z = w_obj*numpy.sqrt(1 + ((wave_length*z*1e-9)/(numpy.pi*w_obj**2))**2)
        N_z = N_0*(1 - numpy.exp(-2*(w_p/w_z)**2))

	# Sscanning System
	# Temporaly : wide field view
        func_r  = numpy.array(map(lambda x : 1.00, r))
        self.source_flux = numpy.array(map(lambda x, y : (2*x)/(numpy.pi*y**2)*func_r, N_z, w_z))



    def set_Detection_path(self) :

        wave_length = self.psf_wavelength*1e-9

        # Magnification
        Mag = 1.0

        # (2) objective lens
        f_obj = self.objective_focal_length

        # (3) tube lens
        f_t = self.tubelens_focal_length1

        # Magnification : Obj to tube lens
        Mag = (f_t/f_obj)*Mag

        # (4) scan lens
        f_s = self.scanlens_focal_length

        # (5) pinhole lens in front of detector
        f_p = self.pinholelens_focal_length

        # Magnification : scan to pinhole lens
        Mag = (f_p/f_s)*Mag

        # set image scaling factor
        voxel_radius = self.spatiocyte_VoxelRadius

        view = self.detector_pixel_length/(2.0*voxel_radius)
        zoom = self.detector_zoom

        self.image_scaling = view/(Mag*zoom)


        # Detector PSF
        self.set_PSF_detector(Mag)




class ConfocalVisualizer(EPIFMVisualizer) :

	'''
	Confocal Visualization class of e-cell simulator
	'''

	def __init__(self, configs=ConfocalConfigs()) :

		assert isinstance(configs, ConfocalConfigs)
		self.configs = configs

		"""
		Check and create the folder for image file.
		"""
		if not os.path.exists(self.configs.movie_image_file_dir):
		    os.makedirs(self.configs.movie_image_file_dir)
		#else:
		#    for file in os.listdir(self.configs.movie_image_file_dir):
		#	os.remove(os.path.join(self.configs.movie_image_file_dir, file))

                """
                set Secondary Image Size and Boundary
                """
                self.img_width  = int(self.configs.detector_image_size[0])
                self.img_height = int(self.configs.detector_image_size[1])

                if self.img_width > IMAGE_SIZE_LIMIT or self.img_height > IMAGE_SIZE_LIMIT :
                        raise VisualizerErrror('Image size is bigger than the limit size')


                """
                Optical Path
                """
		self.configs.set_Optical_path()



        def output_frames(self, num_div=1) :

                # define observational image plane in nm-scale
                voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9

                Nz = int(self.configs.spatiocyte_lengths[0][2]*voxel_size)
                Nx = int(self.configs.spatiocyte_lengths[0][1]*voxel_size)
                Ny = int(self.configs.spatiocyte_lengths[0][0]*voxel_size)

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
                delta_frame = int(round(frame_interval/delta_data))

                # create frame data composed by frame element data
		count = int(round(time/frame_interval))
		count0 = count

                # set number of processors
                max_runs = 200
                #n_cores = multiprocessing.cpu_count()

		while (time < end) :

                    # set image file name
                    image_file_name = os.path.join(self.configs.movie_image_file_dir, \
                                                self.configs.movie_image_file_name_format % (count))

                    print 'time : ', time, ' sec (', count, ')'

                    # define cell
                    #cell = numpy.zeros(shape=(Nz, Ny))
                    mp_arr = multiprocessing.Array(ctypes.c_double, Nz*Ny)
                    np_arr = numpy.frombuffer(mp_arr.get_obj())
                    cell = np_arr.reshape((Nz,Ny))

                    count_start = (count - count0)*delta_frame
                    count_end = (count - count0 + 1)*delta_frame

		    frame_data = self.configs.spatiocyte_data[count_start:count_end]

		    # loop for frame data
                    for i in range(len(frame_data)) :

			# i-th data in a frame
			i_time = frame_data[i][0]
			data   = frame_data[i][1]
			total  = len(data)

			Norm = 1.0#numpy.exp(-(i_time - time)/exposure_time)

                        # loop for particles (multiprocessing)
                        jobs = []

                        for j in range(total) :
                            proc = multiprocessing.Process(target=self.get_molecule_plane, args=(cell, i_time, data[j], p_0, p_b))
                            jobs.append(proc)

                        run = 0
                        while (run < total) :

                            for j in range(max_runs) :

                                if (run + j < total) :
                                    jobs[run+j].start()
                                    sleep(0.1)

                            for j in range(max_runs) :

                                if (run + j < total) :
                                    jobs[run+j].join()

                            run += max_runs


                    if (numpy.amax(cell) > 0) :

			# Output image : Assuming CCD camera output
			camera = self.detector_output_CCD(cell)
			camera.astype('uint%d' % (self.configs.detector_ADC_bit))

			# save image to file
			toimage(camera, low=numpy.amin(camera), high=numpy.amax(camera), mode='I').save(image_file_name)

			# PMT Detector output : intensity (ADC count)
			#intensity = self.detector_output_PMT(cell)

		    time  += frame_interval
		    count += 1



	def detector_output_CCD(self, cell) :

		# Detector Output
		voxel_radius = self.configs.spatiocyte_VoxelRadius
                voxel_size = (2.0*voxel_radius)/1e-9
		re_scaling = (self.configs.image_scaling/self.configs.detector_pixel_length)*(2.0*self.configs.pinholelens_radius)

		Nw_pixel = self.img_width
		Nh_pixel = self.img_height

		Np = int(re_scaling*voxel_size)

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
                cell_pixel = numpy.zeros(shape=(Nw_cell/Np, Nh_cell/Np))

		print cell_pixel.shape

		for i in range(Nh_cell/Np) :
		    for j in range(Nw_cell/Np) :

			# get signal
			signal = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])

			# get background
			background = 0

			# get noise
			noise = self.get_noise(signal + background)

			# get signal + noise
			M  = self.configs.detector_emgain
			#PE = numpy.random.poisson(M*(signal + background + noise), None)
			PE = numpy.random.normal(M*(signal + background), noise, None)

			ADC = self.A2D_converter(PE)
			cell_pixel[i][j] = ADC


                # flat image in pixel-scale
		signal = 0
		background = 0
		noise  = self.get_noise(signal + background)

                #PE = numpy.random.poisson(M*(signal + background + noise), Nw_pixel*Nh_pixel)
                PE = numpy.random.normal(M*(signal + background), noise, Nw_pixel*Nh_pixel)
		camera = numpy.array(map(lambda x : self.A2D_converter(x), PE))
		camera_pixel = camera.reshape([Nw_pixel, Nh_pixel])

		#camera_pixel = numpy.zeros(shape=(Nw_pixel, Nh_pixel))
                #ADC0 = self.configs.detector_ADC_offset
                #camera_pixel.fill(ADC0)

		w_cam_from = int(w_cam_from/Np)
		w_cam_to   = int(w_cam_to/Np)
		h_cam_from = int(h_cam_from/Np)
		h_cam_to   = int(h_cam_to/Np)

		w_cel_from = int(w_cel_from/Np)
		w_cel_to   = int(w_cel_to/Np)
		h_cel_from = int(h_cel_from/Np)
		h_cel_to   = int(h_cel_to/Np)

		ddw = (w_cam_to - w_cam_from) - (w_cel_to - w_cel_from)
		ddh = (h_cam_to - h_cam_from) - (h_cel_to - h_cel_from)

		if   (ddw > 0) : w_cam_to = w_cam_to - ddw
		elif (ddw < 0) : w_cel_to = w_cel_to - ddw

                if   (ddh > 0) : h_cam_to = h_cam_to - ddh
                elif (ddh < 0) : h_cel_to = h_cel_to - ddh


		camera_pixel[w_cam_from:w_cam_to, h_cam_from:h_cam_to] = cell_pixel[w_cel_from:w_cel_to, h_cel_from:h_cel_to]

		return camera_pixel

#		# pinhole in pixel-scale
#                N_ph = (self.configs.image_scaling*voxel_size)
#                PH_radius = 4.0*N_ph/Np
#
#                w = numpy.linspace(0, Nw_pixel-1, Nw_pixel)
#                h = numpy.linspace(0, Nh_pixel-1, Nh_pixel)
#		W, H = numpy.meshgrid(w, h)
#
#		y0 = Nw_pixel/2.0
#		z0 = Nh_pixel/2.0
#
#		func = lambda y, z : abs(PH_radius**2 - ((y - y0)**2 + (z - z0)**2))
#		pinhole_pixel = numpy.array(map(lambda z : map(lambda y : 4e+5 if func(y, z) < 3 else 0.00, w), h))
#
#		# add pinhole to camera
#		camera_pixel += pinhole_pixel
#
#		return camera_pixel



	def detector_output_PMT(self, cell) :

                voxel_size = (2.0*self.configs.spatiocyte_VoxelRadius)/1e-9

                Nw_pixel = 1#self.img_width
                Nh_pixel = 1#self.img_height

                # pinhole size
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

		# get signal through the pinhole
		signal = numpy.sum(plane)

                # get background
                background = 0

                # get noise
                noise = self.get_noise(signal + background)

                # get signal + noise
                M  = self.configs.detector_emgain
                #PE = numpy.random.poisson(M*(signal + background + noise), None)
                PE = numpy.random.normal(M*(signal + background), noise, None)

                # convert signals in photoelectron to count
                ADC = self.A2D_converter(PE)

		return ADC


