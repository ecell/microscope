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


class FCSConfigs(EPIFMConfigs) :

    '''
    FCS configuration

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

#        # (2) beam expander
#        f_1 = self.expander_focal_length1
#        f_2 = self.expander_focal_length2
#
#        w_p = self.expander_pinhole_radius
#
#        w_BE = (f_2/f_1)*w_source
#
#        # (3) scan and tube lens
#        f_s = self.scanlens_focal_length
#        f_t = self.tubelens_focal_length1
#
#        w_tube = (f_t/f_s)*w_BE
#
#        # (4) objective
#        f_obj = self.objective_focal_length
#
#        # Rayleigh range
#        z_R = numpy.pi*w_tube**2/wave_length
#
#        # object distance to maximize image distance
#        s_obj = f_obj + z_R
#        w_obj = w_tube/numpy.sqrt((1 - s_obj/f_obj)**2 + (z_R/f_obj)**2)

        # Beam Flux [photons/(m^2 sec)]
	w_source = 200e-9
        w_z = w_obj*numpy.sqrt(1 + ((wave_length*z*1e-9)/(numpy.pi*w_source**2))**2)
        #N_z = N_0*(1 - numpy.exp(-2*(w_p/w_z)**2))

        self.source_flux = numpy.array(map(lambda x : (2*N0)/(numpy.pi*x**2)*numpy.exp(-2*(r*1e-9/x)**2), w_z))

	# Temporaly : wide field view
        #func_r  = numpy.array(map(lambda x : 1.00, r))
        #self.source_flux = numpy.array(map(lambda x, y : (2*x)/(numpy.pi*y**2)*func_r, N_z, w_z))



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

        view = self.pinholelens_radius/voxel_radius
        zoom = self.detector_zoom

        self.pinhole_scaling = view/(Mag*zoom)


        # Detector PSF
        self.set_PSF_detector(Mag)




class FCSVisualizer(EPIFMVisualizer) :

	'''
	FCS Visualization class of e-cell simulator
	'''

	def __init__(self, configs=FCSConfigs()) :

		assert isinstance(configs, FCSConfigs)
		self.configs = configs

		"""
		Check and create the folder for image file.
		"""
		if not os.path.exists(self.configs.movie_image_file_dir):
		    os.makedirs(self.configs.movie_image_file_dir)
		#else:
		#    for file in os.listdir(self.configs.movie_image_file_dir):
		#	os.remove(os.path.join(self.configs.movie_image_file_dir, file))

                if not os.path.exists(self.configs.output_file_dir):
                    os.makedirs(self.configs.output_file_dir)

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



	def get_noise(self, signal) :

		# detector noise in current unit (Ampare)
                NA = self.configs.detector_readout
                Id = self.configs.detector_dark_current
                Fn = numpy.sqrt(self.configs.detector_excess)
                B  = self.configs.detector_bandwidth
                M  = self.configs.detector_gain
                e  = self.configs.electron_charge

                # Noise calculation defined in Hamamatsu PMT technical guide
                sigma2 = 2*e*B*Fn**2*M**2*(signal + Id)+ (NA)**2
                noise  = numpy.sqrt(sigma2)
        
                return noise



        def get_signal(self, time, pid, s_index, p_i, p_b, p_0) :

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

                if (d_s < len(d)) :
		    source_depth = d_s
                else :
		    source_depth = d[-1]

		# beam lateral position
		rr = numpy.sqrt((y_i-y_0)**2 + (z_i-z_0)**2)

		if (rr < len(r)) :
		    source_radius = rr
		else :
		    source_radius = r[-1]

                # normalization
                radius = 1e-9
		unit_area = radius**2
		norm = unit_area

		# get illumination PSF
		source_psf = norm*self.configs.source_flux[int(source_depth)][int(source_radius)]
		source_max = norm*self.configs.source_flux[0][0]

                # signal conversion : Output Intensity = Physics * PSF (Beam)
		#Intensity = self.get_intensity(time, pid, source_psf, source_max)
		Ratio = self.effects.conversion_ratio
		Intensity = Ratio * source_psf

		# fluorophore axial position
		if (d_f < len(d)) :
		    fluo_depth = d_f
		else :
		    fluo_depth = d[-1]

		# coordinate transformation : polar --> cartisian
		theta = numpy.linspace(0, 180, 181)

                z = numpy.linspace(0, +r[-1], len(r))
                y = numpy.linspace(-r[-1], +r[-1], 2*len(r)-1)

		psf_t = numpy.array(map(lambda x : 1.00, theta))
		psf_r = self.configs.fluorophore_psf[int(fluo_depth)]

		psf_polar = numpy.array(map(lambda x : psf_t*x, psf_r))

                # get fluorophore PSF
		fluo_psf  = numpy.array(self.polar2cartesian(r, theta, psf_polar, z, y))

                # signal conversion : Output PSF = Intensity * PSF(Fluorophore)
		signal = Intensity * fluo_psf


                return signal



        def output_frames(self, num_div=1) :

                # define observational image plane in nm-scale
                voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9

                Nz = int(self.configs.spatiocyte_lengths[0][2]*voxel_size)
                Ny = int(self.configs.spatiocyte_lengths[0][1]*voxel_size)
                Nx = int(self.configs.spatiocyte_lengths[0][0]*voxel_size)

                # focal point
                p_0 = numpy.array([Nx, Ny, Nz])*self.configs.detector_focal_point

                # beam position : Assuming beam position = focal point for temporary
                p_b = copy.copy(p_0)

		# frame interval
		frame_interval = 1.0/self.configs.detector_fps
		#frame_interval = self.configs.detector_frame_interval
		exposure_time  = self.configs.detector_exposure_time

		time = self.configs.detector_start_time
		end  = self.configs.detector_end_time

		# data-time interval
                t0 = self.configs.spatiocyte_data[0][0]
                t1 = self.configs.spatiocyte_data[1][0]

		delta_data  = t1 - t0
                delta_time = int(round(exposure_time/delta_data))
                #delta_time = int(round(exposure_time/frame_interval))

                # create frame data composed by frame element data
		count  = int(round(time/exposure_time))
		count0 = count

		# initialize Physical effects
		#length0 = len(self.configs.spatiocyte_data[0][1])
		#self.effects.set_states(t0, length0)

		# set number of processors
		max_runs = 100
		#max_runs = multiprocessing.cpu_count()

		# set A/D Converter
		self.set_ADConverter()

		while (time < end) :

                    # set image file name
                    self.image_file_name = os.path.join(self.configs.movie_image_file_dir, \
						self.configs.movie_image_file_name_format % (count))

                    print 'time : ', time, ' sec (', count, ')'

		    # define cell
                    #cell = numpy.zeros(shape=(Nz, Ny))
		    mp_arr = multiprocessing.Array(ctypes.c_double, Nz*Ny)
		    np_arr = numpy.frombuffer(mp_arr.get_obj())
		    cell = np_arr.reshape((Nz, Ny))

		    count_start = (count - count0)*delta_time
		    count_end   = (count - count0 + 1)*delta_time

		    frame_data = self.configs.spatiocyte_data[count_start:count_end]

		    # loop for frame data
		    for i in range(len(frame_data)) :

			# i-th data in a frame
			i_time = frame_data[i][0]
			data   = frame_data[i][1]
			total  = len(data)

                        print '\t', '%02d-th frame : ' % (i), i_time, ' sec'

			# loop for particles (multiprocessing)
			jobs = []

			for j in range(total) :
			    proc = multiprocessing.Process(target=self.get_molecule_plane, args=(cell, i_time-time, data[j], j, p_0, p_b))
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

			camera = self.detector_output(p0, cell)
			camera.astype('uint%d' % (self.configs.detector_ADC_bit))

			# save image to file
			#toimage(camera, low=numpy.amin(camera), high=numpy.amax(camera), mode='I').save(self.image_file_name)
			toimage(camera, cmin=1600, cmax=6000).save(self.image_file_name)
			#scipy.misc.imsave(self.image_file_name, camera)


		    time  += exposure_time
		    count += 1



	def detector_output(self, p0, cell) :

		# get focal position
		x0, y0, z0 = p0

                # Pinhole randius
		voxel_radius = self.configs.spatiocyte_VoxelRadius
		pinhole_radius = pinhole_scaling*voxel_radius/1e-9

                # image in nm-scale
		Nw_cell = len(cell)
		Nh_cell = len(cell[0])

		# get photon flux (Photons/sec)
		photon_flux= 0

		for i in range(Nw_cell) :
		    for j in range(Nh_cell) :

			# get photons
			if ((i - y0)**2 + (j - z0)**2 < pinhole_radius**2) :
			    photon_flux += cell[i][j]


		# get Quantum Efficiency
		index = int(self.configs.psf_wavelength) - int(self.configs.wave_length[0])
		QE = self.configs.detector_qeff[index]

		# get electron charge
		e = 1.62e-19 # coulomb
		#e = self.configs.electron_charge

                # get current at cathode (photoelectric current)
		cathode_current = e*QE*photon_flux

                # get detector noise (photoelectric current)
                noise = self.get_noise(cathode_current + background)

                # get Gain
                M  = self.configs.detector_emgain

                # get current at anode (photoelectric current)
		anode_current = numpy.random.normal(M*(cathode_current + background), noise, None)

		# Gated integrator
		charge = anode_current * interval

		# A/D converter : ADC counts
		ADC = self.get_ADC_value(charge)

		return ADC




	def set_ADConverter(self) :

	    # set ADC offset
	    ADC_offset = self.configs.detector_ADC_offset

	    # set ADC gain
            ADC_bit  = self.configs.detector_ADC_bit
            ADC_satQ = self.configs.detector_ADC_satQ

            ADC_gain = (ADC_satQ - 0)/(2**ADC_bit - ADC_offset)

	    print 'A/D Converter :%d-bit' % (int(ADC_bit))
	    print '\toffset =', ADC_offset
	    print '\tgain   =', ADC_gain

            self.configs.detector_ADC_gain   = ADC_gain



	def get_ADC_value(self, charge) :

	    # check non-linearity
	    Q_max = self.configs.detector_ADC_satQ

	    if (photo_electron > Q_max) :
		photo_electron = Q_max

            # convert photoelectron to ADC counts (Grayscale)
            k = self.configs.detector_ADC_gain
            ADC0 = self.configs.detector_ADC_offset
            ADC_max = 2**self.configs.detector_ADC_bit - 1

            ADC = charge/k + ADC0

	    if (ADC > ADC_max) :
		ADC = ADC_max

	    if (ADC < 0) :
		ADC = 0

	    return int(ADC)



