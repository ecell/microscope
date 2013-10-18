"""

linescan_confm_handler.py

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
import ctypes
import multiprocessing

import scipy
import numpy

import parameter_configs
from effects_handler import PhysicalEffects
from epifm_handler import VisualizerError, EPIFMConfigs, EPIFMVisualizer

from scipy.special import j0
from scipy.misc    import toimage


class LineScanConfocalConfigs(EPIFMConfigs) :

    '''
    Line-scanning Confocal configuration

	Line-like Gaussian Profile
	    +
	Line-scanning
	    +
	Slit
	    +
	Detector : EMCCD
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



    def set_Slits(self, size = None) :

        print '--- Slit :'

        self._set_data('slit_size', size)

        print '\tSize = ', self.slit_size, 'm'



    def set_Illumination_path(self) :

        r = self.radial
        d = numpy.linspace(0, 20000, 20001)
        #d = self.depth

        # (plank const) * (speed of light) [joules meter]
        hc = 2.00e-25

	# Illumination
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

	# Rayleigh range
        z_R = numpy.pi*w_0**2/wave_length

        # Beam Flux [photons/(m^2 sec)]
        w_z = w_0*numpy.sqrt(1 + ((wave_length*d*1e-9)/(numpy.pi*w_0**2))**2)

	# photon flux density [photon/(sec m^2)]
        self.source_flux = numpy.array(map(lambda x : 2*N_0/(numpy.pi*x**2)*numpy.exp(-2*(r*1e-9/x)**2), w_z))

	print 'Photon Flux Density (Max) :', numpy.amax(self.source_flux)



    def set_Detection_path(self) :

        wave_length = self.psf_wavelength*1e-9

	# Magnification
	Mag = self.image_magnification

	# set image scaling factor
        voxel_radius = self.spatiocyte_VoxelRadius

	# set zoom
	zoom = self.detector_zoom

	# set slit pixel length
	pixel_length = self.slit_size/(Mag*zoom)

	self.image_resolution = pixel_length
	self.image_scaling = pixel_length/(2.0*voxel_radius)

	print 'Resolution :', self.image_resolution, 'm'
	print 'Scaling :', self.image_scaling

        # Detector PSF
        self.set_PSF_detector()



class LineScanConfocalVisualizer(EPIFMVisualizer) :

	'''
	Confocal Visualization class of e-cell simulator
	'''

	def __init__(self, configs=LineScanConfocalConfigs(), effects=PhysicalEffects()) :

		assert isinstance(configs, LineScanConfocalConfigs)
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
                Optical Path
                """
		self.configs.set_Optical_path()



        def get_signal(self, time, pid, s_index, p_i, p_b, p_0, norm) :

		# set focal point
		x_0, y_0, z_0 = p_0

                # set source center
                x_b, y_b, z_b = p_b

		# set particle position
                x_i, y_i, z_i = p_i

		#
		r = self.configs.radial
		d = self.configs.depth

		# beam axial position
		d_s = abs(x_i - x_b)

                if (d_s < 20000) :
		    source_depth = d_s
                else :
		    source_depth = 19999

		# beam horizontal position (y-direction)
		hh = numpy.sqrt((y_i-y_b)**2)

		if (hh < len(r)) :
		    source_horizon = hh
		else :
		    source_horizon = r[-1]

		# get illumination PSF
		source_psf = self.configs.source_flux[int(source_depth)][int(source_horizon)]
		#source_max = norm*self.configs.source_flux[0][0]

                # signal conversion : 
		#Intensity = self.get_intensity(time, pid, source_psf, source_max)
		Ratio = self.effects.conversion_ratio

		# fluorophore axial position
		d_f = abs(x_i - x_b)

		if (d_f < len(d)) :
		    fluo_depth = d_f
		else :
		    fluo_depth = d[-1]

		# get fluorophore PSF
		fluo_psf = self.fluo_psf[int(fluo_depth)]

                # signal conversion : Output PSF = PSF(source) * Ratio * PSF(Fluorophore)
		signal = norm * source_psf * Ratio * fluo_psf

                return signal



	def get_molecule_plane(self, cell, time, data, pid, p_b, p_0) :

		voxel_size = (2.0*self.configs.spatiocyte_VoxelRadius)/1e-9

		# get beam position
		x_b, y_b, z_b = p_b

                # cutoff randius
		slit_radius = int(self.configs.image_scaling*voxel_size/2)
		cut_off = 3*slit_radius

		# particles coordinate, species and lattice IDs
                c_id, s_id, l_id = data

		sid_array = numpy.array(self.configs.spatiocyte_species_id)
		s_index = (numpy.abs(sid_array - int(s_id))).argmin()

		if self.configs.spatiocyte_observables[s_index] is True :

		    # Normalization
		    unit_time = 1.0
		    unit_area = (1e-9)**2
		    norm = (unit_area*unit_time)/(4.0*numpy.pi)

		    # particles coordinate in real(nm) scale
                    #pos = self.get_coordinate(c_id)
                    #p_i = numpy.array(pos)*voxel_size
                    p_i = self.get_coordinate(c_id)
		    x_i, y_i, z_i = p_i

		    if (abs(y_i - y_b) < cut_off) :

                        #print pid, s_id, p_i
                        # get signal matrix
                        signal = self.get_signal(time, pid, s_index, p_i, p_b, p_0, norm)

                        # add signal matrix to image plane
                        self.overwrite_signal(cell, signal, p_i)



	def output_frames(self, num_div=1):

	    # set Fluorophores PSF
	    self.set_fluo_psf()

            start = self.configs.spatiocyte_start_time
            end = self.configs.spatiocyte_end_time

            exposure_time = self.configs.detector_exposure_time
            num_timesteps = int(math.ceil((end - start) / exposure_time))

	    index0  = int(round(start/exposure_time))

            envname = 'ECELL_MICROSCOPE_SINGLE_PROCESS'

            if envname in os.environ and os.environ[envname]:
                self.output_frames_each_process(index0, num_timesteps)
            else:
                num_processes = multiprocessing.cpu_count()
                n, m = divmod(num_timesteps, num_processes)
                # when 10 tasks is distributed to 4 processes,
                # number of tasks of each process must be [3, 3, 2, 2]
                chunks = [n + 1 if i < m else n for i in range(num_processes)]

                processes = []
                start_index = index0

                for chunk in chunks:
                    stop_index = start_index + chunk
                    process = multiprocessing.Process(
                        target=self.output_frames_each_process,
                        args=(start_index, stop_index))
                    process.start()
                    processes.append(process)
                    start_index = stop_index

                for process in processes:
                    process.join()



        def output_frames_each_process(self, start_count, stop_count):

		voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9

		# image dimenssion in pixel-scale
		Nw_pixel = self.configs.detector_image_size[0]
		Nh_pixel = self.configs.detector_image_size[1]

		# cells dimenssion in nm-scale
                Nz = int(self.configs.spatiocyte_lengths[2] * voxel_size)
                Ny = int(self.configs.spatiocyte_lengths[1] * voxel_size)
                Nx = int(self.configs.spatiocyte_lengths[0] * voxel_size)

		# pixel length : nm/pixel
		Np = int(self.configs.image_scaling*voxel_size)

		# cells dimenssion in pixel-scale
		Ny_pixel = Ny/Np
		Nz_pixel = Nz/Np

                # focal point
                p_0 = numpy.array([Nx, Ny, Nz])*self.configs.detector_focal_point

                # beam position : 
                beam_center = numpy.array(self.configs.detector_focal_point)

                # set boundary condition
                if (self.configs.spatiocyte_bc_switch == True) :

                    bc = numpy.zeros(shape=(Nz, Ny))
                    bc = self.set_boundary_plane(bc, p_b, p_0)

		# exposure time on cell
		R = float(Ny_pixel) / float(Nw_pixel)

                exposure_time = self.configs.detector_exposure_time
                contact_time  = R * exposure_time
		non_contact_time = (exposure_time - contact_time)/2

                spatiocyte_start_time = self.configs.spatiocyte_start_time
                time = exposure_time * start_count
                end  = exposure_time * stop_count

                # data-time interval
		data_interval = self.configs.spatiocyte_interval

		# time/count
                delta_time = int(round(exposure_time / data_interval))

                # create frame data composed by frame element data
                count  = start_count
                count0 = int(round(spatiocyte_start_time / exposure_time))

                # initialize Physical effects
                #length0 = len(self.configs.spatiocyte_data[0][1])
                #self.effects.set_states(t0, length0)

                while (time < end) :

                    # set image file name
                    image_file_name = os.path.join(self.configs.image_file_dir,
                                        self.configs.image_file_name_format % (count))

                    print 'time : ', time, ' sec (', count, ')'

                    # define cell
                    cell = numpy.zeros(shape=(Nz, Ny))

                    count_start = (count - count0)*delta_time
                    count_end   = (count - count0 + 1)*delta_time

                    frame_data = self.configs.spatiocyte_data[count_start:count_end]

		    if (len(frame_data) > 0) :

		        # beam position : initial
		        p_b = numpy.array([Nx, Ny, Nz])*beam_center

		        # line-scanning sequences
		        scan_time = 0

		        while (scan_time < contact_time) :

			    # define : cell at scanned region
			    scan_cell = numpy.zeros(shape=(Nz, Ny))

			    # beam position : i-th frame
			    if (count%2 == 0) :
				# scan from left-to-right
				beam_center[1] = scan_time/contact_time
			    else :
				# scan from right-to-left
				beam_center[1] = 1 - scan_time/contact_time

			    p_b = numpy.array([Nx, Ny, Nz])*beam_center
			    x_b, y_b, z_b = p_b

                            # loop for frame data
			    i_time, i_data = frame_data[0]

			    diff = abs(i_time - (scan_time + time + non_contact_time))
			    data = i_data

                            for i, (i_time, i_data) in enumerate(frame_data) :
                                #print '\t', '%02d-th frame : ' % (i), i_time, ' sec'
			        i_diff = abs(i_time - (scan_time + time))

			        if (i_diff < diff) :
				    diff = i_diff
				    data = i_data

                            # loop for particles
                            for j, j_data in enumerate(data) :
                                self.get_molecule_plane(scan_cell, i_time, j_data, j, p_b, p_0)

			    # overwrite the scanned region to cell
			    r_s = int(self.configs.image_scaling*voxel_size/2)

			    if (y_b-r_s < 0) : y_from = int(y_b)
			    else : y_from = int(y_b - r_s)

			    if (y_b+r_s > Ny) : y_to = int(y_b)
			    else : y_to = int(y_b + r_s)

			    cell[:, y_from:y_to] += scan_cell[:, y_from:y_to]*(contact_time/Ny_pixel)

			    scan_time += contact_time/Ny_pixel


		    if (numpy.amax(cell) > 0) :

			if (self.configs.spatiocyte_bc_switch == True) :
			    camera = self.detector_output(cell, bc)
			else : camera = self.detector_output(cell)

			camera.astype('uint%d' % (self.configs.ADConverter_bit))

			# save image to file
			toimage(camera, low=numpy.amin(camera), high=numpy.amax(camera), mode='I').save(image_file_name)

		    time  += exposure_time
		    count += 1



	def detector_output(self, cell, bc=None) :

		# Detector Output
		voxel_radius = self.configs.spatiocyte_VoxelRadius
                voxel_size = (2.0*voxel_radius)/1e-9

		Nw_pixel = self.configs.detector_image_size[0]
		Nh_pixel = self.configs.detector_image_size[1]

		# set camera's image scaling
		Mag  = self.configs.image_magnification
		zoom = self.configs.detector_zoom
		view = self.configs.detector_pixel_length/(2.0*voxel_radius)

		camera_image_scaling = view/(Mag*zoom)

		Np = int(camera_image_scaling*voxel_size)

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

		# set seed for random number
		numpy.random.seed()

		# Photon distribution
		for i in range(Nw_cell/Np) :
		    for j in range(Nh_cell/Np) :

			# get photons
			photons = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])

			if (photons > 0) :

			    # get crosstalk
			    if (self.effects.detector_crosstalk_switch == True) :

				width = self.effects.detector_crosstalk_width

				n_i = numpy.random.normal(0, width, photons)
				n_j = numpy.random.normal(0, width, photons)

				#i_bins = int(numpy.amax(n_i) - numpy.amin(n_i))
				#j_bins = int(numpy.amax(n_j) - numpy.amin(n_j))

				smeared_photons, edge_i, edge_j = numpy.histogram2d(n_i, n_j, bins=(24, 24), range=[[-12,12],[-12,12]])

				# smeared photon distributions
				cell_pixel = self.overwrite_smeared(cell_pixel, smeared_photons, i, j)

			    else :

				cell_pixel[i][j] = photons

		index = int(self.configs.psf_wavelength) - int(self.configs.wave_length[0])
		QE = self.configs.detector_qeff[index]

		# Photoelectron and ADC count distribution
                for i in range(Nw_cell/Np) :
                    for j in range(Nh_cell/Np) :

			# pixel position
			pixel = (i, j)

			# Detector : Quantum Efficiency
			index = int(self.configs.psf_wavelength) - int(self.configs.wave_length[0])
			QE = self.configs.detector_qeff[index]

                        # get signal (photoelectrons)
			signal = QE*cell_pixel[i][j]

                        # get constant background (photoelectrons)
                        if (self.effects.background_switch == True) :

                            mean = self.effects.background_mean
                            background = QE*numpy.random.poisson(mean, None)

                        else : background = 0

                        # get detector noise (photoelectrons)
                        noise = self.get_noise(signal + background)

                        # get EM Gain
                        M = self.configs.detector_emgain

                        # get signal + background (photoelectrons)
			PE = numpy.random.normal(M*(signal + background), noise, None)

			# A/D converter : Photoelectrons --> ADC counts
			ADC = self.get_ADC_value(pixel, PE)

			if (self.configs.spatiocyte_bc_switch == True) :
			    cell_pixel[i][j] = ADC*bc[i][j]
			else :
			    cell_pixel[i][j] = ADC


                # Background (No photon signal)
                camera_pixel = numpy.zeros([Nw_pixel, Nh_pixel])

                for i in range(Nw_pixel) :
                    for j in range(Nh_pixel) :

			# set pixel position
			pixel = (i, j)

			signal, background = 0, 0
			noise = self.get_noise(signal + background)

			PE = numpy.random.normal(M*(signal + background), noise, None)
			#PE = M*(signal + background)
			ADC = self.get_ADC_value(pixel, PE)
			camera_pixel[i][j] = ADC
			#camera_pixel[i][j] = PE


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



