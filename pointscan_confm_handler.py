"""

linescan_confm_handler.py

"""

import sys
import os
import copy
import math
import multiprocessing

import numpy

import parameter_configs
from effects_handler import PhysicalEffects
from fcs_handler import FCSConfigs, FCSVisualizer

from scipy.misc import toimage


class PointScanConfocalConfigs(FCSConfigs) :

    '''
    Point-scanning Confocal configuration

	Point-like Gaussian Profile
	    +
	Point-scanning
	    +
	Pinhole
	    +
	Detector : PMT
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



    def set_Pinhole(self, radius = None) :

        print '--- Pinhole :'

        self._set_data('pinhole_radius', radius)

        print '\tRadius = ', self.pinhole_radius, 'm'



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

	# set pinhole pixel length
	pixel_length = (2.0*self.pinhole_radius)/(Mag*zoom)

	self.image_resolution = pixel_length
	self.image_scaling = pixel_length/(2.0*voxel_radius)

	print 'Magnification : x %d' % (Mag)
	print 'Resolution :', self.image_resolution, 'm'
	print 'Scaling :', self.image_scaling

        # Detector PSF
        self.set_PSF_detector()




class PointScanConfocalVisualizer(FCSVisualizer) :

	'''
	Confocal Visualization class of e-cell simulator
	'''

	def __init__(self, configs=PointScanConfocalConfigs(), effects=PhysicalEffects()) :

		assert isinstance(configs, PointScanConfocalConfigs)
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

        def get_signal(self, time, pid, s_index, p_i, p_b, p_0, norm):

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

            if d_s < 20000:
                source_depth = d_s
            else:
                source_depth = 19999

            # beam horizontal position (y-direction)
            hh = numpy.sqrt((y_i-y_b)**2)

            if hh < len(r):
                source_horizon = hh
            else:
                source_horizon = r[-1]

            # get illumination PSF
            source_psf = self.configs.source_flux[int(source_depth)][int(source_horizon)]
            #source_max = norm*self.configs.source_flux[0][0]

            # signal conversion :
            #Intensity = self.get_intensity(time, pid, source_psf, source_max)
            Ratio = self.effects.conversion_ratio

            # fluorophore axial position
            d_f = abs(x_i - x_b)

            if d_f < len(d):
                fluo_depth = d_f
            else:
                fluo_depth = d[-1]

            # get fluorophore PSF
            fluo_psf = self.fluo_psf[int(fluo_depth)]

            # signal conversion : Output PSF = PSF(source) * Ratio * PSF(Fluorophore)
            signal = norm * source_psf * Ratio * fluo_psf

            return signal

        def get_molecule_plane(self, cell, time, data, pid, p_b, p_0, offset):
            voxel_size = (2.0*self.configs.spatiocyte_VoxelRadius)/1e-9

            # get beam position
            x_b, y_b, z_b = p_b

            # cutoff randius
            pinhole_radius = int(self.configs.image_scaling*voxel_size/2)
            cut_off = 3*pinhole_radius

            # particles coordinate, species and lattice IDs
            c_id, s_id, l_id = data

            sid_array = numpy.array(self.configs.spatiocyte_species_id)
            s_index = (numpy.abs(sid_array - int(s_id))).argmin()

            if self.configs.spatiocyte_observables[s_index] is True:

                # Normalization
                unit_time = 1.0
                unit_area = (1e-9)**2
                norm = (unit_area*unit_time)/(4.0*numpy.pi)

                # particles coordinate in nm-scale
                #pos = self.get_coordinate(c_id)
                #p_i = numpy.array(pos)*voxel_size
                p_i = self.get_coordinate(c_id)
                x_i, y_i, z_i = p_i

                if numpy.sqrt((y_i - y_b)**2 + (z_i - z_b)**2) < cut_off:

                    #print pid, s_id, p_i
                    # get signal matrix
                    signal = self.get_signal(time, pid, s_index, p_i, p_b, p_0, norm)

                    z_from, y_from = offset
                    z_size, y_size = cell.shape
                    z_to = z_from + y_size
                    y_to = y_from + y_size
                    r = len(self.configs.radial)

                    self.overwrite_signal(cell, signal, p_i, r, z_from, y_from, z_to, y_to)

        def overwrite_signal(self, cell, signal, p_i, r, z_from, y_from, z_to, y_to):
            z_size = z_to - z_from
            y_size = y_to - y_from
            x_i, y_i, z_i = p_i
            zi_size, yi_size = signal.shape

            z0_from = bounded(z_i - r - z_from, lower_bound=0, upper_bound=z_size)
            y0_from = bounded(y_i - r - y_from, lower_bound=0, upper_bound=y_size)
            z0_to = bounded(z_i + r - z_to, lower_bound=0, upper_bound=z_size)
            y0_to = bounded(y_i + r - y_to, lower_bound=0, upper_bound=y_size)

            zi_from = bounded(z_from - z_i + r, lower_bound=0, upper_bound=zi_size)
            yi_from = bounded(y_from - y_i + r, lower_bound=0, upper_bound=yi_size)
            zi_to = bounded(z_from - z_i + r + z_size, lower_bound=0, upper_bound=zi_size)
            yi_to = bounded(y_from - y_i + r + y_size, lower_bound=0, upper_bound=yi_size)

            cell[z0_from:z0_to, y0_from:y0_to] += signal[zi_from:zi_to, yi_from:yi_to]

        def output_frames(self, num_div=1):

            # set Fluorophores PSF
            self.set_fluo_psf()

            start = self.configs.spatiocyte_start_time
            end = self.configs.spatiocyte_end_time

            exposure_time = self.configs.detector_exposure_time
            num_timesteps = int(math.ceil((end - start) / exposure_time))

            index0 = int(round(start/exposure_time))

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

            # Beam position
            beam_center = numpy.array(self.configs.detector_focal_point)

            # set boundary condition
            if self.configs.spatiocyte_bc_switch is True:
                bc = numpy.zeros(shape=(Nz, Ny))
                bc = self.set_boundary_plane(bc, p_b, p_0)

            # exposure time
            exposure_time = self.configs.detector_exposure_time

            # contacting time with cells vertical-axis
            R_z = float(Nz_pixel) / float(Nh_pixel)

            z_exposure_time = exposure_time
            z_contact_time = R_z*z_exposure_time
            z_nocontact_time = (z_exposure_time - z_contact_time)/2

            # contacting time with cells horizontal-axis
            R_y = float(Ny_pixel) / float(Nw_pixel)

            y_exposure_time = z_exposure_time/Nw_pixel
            y_contact_time = R_y*y_exposure_time
            y_nocontact_time = (y_exposure_time - y_contact_time)/2

            #####
            start_time = self.configs.spatiocyte_start_time

            time = exposure_time * start_count
            end = exposure_time * stop_count

            # data-time interval
            data_interval = self.configs.spatiocyte_interval

            delta_time = int(round(exposure_time / data_interval))

            # create frame data composed by frame element data
            count = start_count
            count0 = int(round(start_time / exposure_time))

            # initialize Physical effects
            #length0 = len(self.configs.spatiocyte_data[0][1])
            #self.effects.set_states(t0, length0)

            while time < end:

                # set image file name
                image_file_name = os.path.join(
                    self.configs.image_file_dir,
                    self.configs.image_file_name_format % (count))

                print 'time : ', time, ' sec (', count, ')'

                # define cell in pixel-scale
                cell = numpy.zeros(shape=(Nz, Ny))

                count_start = (count - count0)*delta_time
                count_end = (count - count0 + 1)*delta_time

                frame_data = self.configs.spatiocyte_data[count_start:count_end]

                if len(frame_data) > 0:

                    # Beam position : initial
                    p_b = numpy.array([Nx, Ny, Nz])*beam_center

                    # contact time : z-direction
                    z_scan_time = 0

                    # no contact time : z-direction
                    non_contact_time = z_nocontact_time

                    # vertical-scanning sequences
                    while z_scan_time < z_contact_time:

                        # Beam position : z-direction
                        beam_center[2] = z_scan_time/z_contact_time

                        # contact time : y-direction
                        y_scan_time = 0

                        # no contact time : y-direction (left-margin)
                        non_contact_time += y_nocontact_time

                        # horizontal-scanning sequences
                        while y_scan_time < y_contact_time:

                            # Beam position : y-direction
                            beam_center[1] = y_scan_time/y_contact_time

                            p_b = numpy.array([Nx, Ny, Nz])*beam_center
                            x_b, y_b, z_b = p_b

                            # loop for frame data
                            i_time, i_data = frame_data[0]

                            scan_time = z_scan_time + y_scan_time + non_contact_time
                            diff = abs(i_time - (scan_time + time))
                            data = i_data

                            for i, (i_time, i_data) in enumerate(frame_data):
                                #print '\t', '%02d-th frame : ' % (i), i_time, ' sec'
                                i_diff = abs(i_time - (scan_time + time))

                                if i_diff < diff:
                                    diff = i_diff
                                    data = i_data

                            # overwrite the scanned region to cell
                            r_p = int(self.configs.image_scaling*voxel_size/2)

                            if y_b-r_p < 0:
                                y_from = int(y_b)
                            else:
                                y_from = int(y_b - r_p)

                            if y_b+r_p >= Ny:
                                y_to = int(y_b)
                            else:
                                y_to = int(y_b + r_p)

                            if z_b-r_p < 0:
                                z_from = int(z_b)
                            else:
                                z_from = int(z_b - r_p)

                            if z_b+r_p >= Nz:
                                z_to = int(z_b)
                            else:
                                z_to = int(z_b + r_p)

                            offset = (z_from, y_from)
                            mask = numpy.zeros(shape=(z_to-z_from, y_to-y_from))

                            zz, yy = numpy.ogrid[z_from-int(z_b):z_to-int(z_b), y_from-int(y_b):y_to-int(y_b)]
                            rr_cut = yy**2 + zz**2 < r_p**2
                            mask[rr_cut] = 1

                            scan_cell = numpy.zeros_like(mask)
                            # loop for particles
                            for j, j_data in enumerate(data):
                                self.get_molecule_plane(scan_cell, i_time, j_data, j, p_b, p_0, offset)

                            cell[z_from:z_to, y_from:y_to] += mask * scan_cell

                            y_scan_time += y_contact_time/Ny_pixel

                        # no contact time : y-direction (right-margin)
                        non_contact_time += y_nocontact_time

                        z_scan_time += z_contact_time/Nz_pixel

                if numpy.amax(cell) > 0:

                    if self.configs.spatiocyte_bc_switch is True:
                        camera = self.detector_output(cell, bc)
                    else:
                        camera = self.detector_output(cell)

                    camera.astype('uint%d' % (self.configs.ADConverter_bit))

                    # save image to file
                    toimage(camera, low=numpy.amin(camera), high=numpy.amax(camera), mode='I').save(image_file_name)

                time += exposure_time
                count += 1



	def detector_output(self, cell, bc=None) :

		# Detector Output
		voxel_radius = self.configs.spatiocyte_VoxelRadius
                voxel_size = (2.0*voxel_radius)/1e-9

		Nw_pixel = self.configs.detector_image_size[0]
		Nh_pixel = self.configs.detector_image_size[1]

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
                cell_pixel = numpy.zeros(shape=(Nw_cell/Np, Nh_cell/Np))

		# set seed for random number
		numpy.random.seed()

		# detector noise in current unit (Ampere)
                NA = self.configs.detector_readout
                Id = self.configs.detector_dark_current
                F  = self.configs.detector_excess
                B  = self.configs.detector_bandwidth
                M  = self.configs.detector_emgain
		e  = self.configs.electron_charge
		D  = Id/e

		# observational time
		T = 1/(2*B)

		# Photoelectron and ADC count distribution
                for i in range(Nw_cell/Np) :
                    for j in range(Nh_cell/Np) :

			# pixel position
			pixel = (i, j)

			# Detector : Quantum Efficiency
			#index = int(self.configs.psf_wavelength) - int(self.configs.wave_length[0])
			#QE = self.configs.detector_qeff[index]
			QE = 0.3

			# get photons
			photon_flux = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])

			# get background
			background = 0

#			if (self.configs.detector_mode == 'Current') :
#
#                	    # get current at cathode (photoelectric current)
#			    cathode_current = e*(QE*photon_flux)
#
#			    # get detector noise (photoelectric current)
#			    noise_current = numpy.random.normal(0, Id, None)
#
#                	    # get anode current (Gain x cathode current)
#			    G = self.configs.detector_emgain
#			    anode_current = G*cathode_current
#
#			    # get total photoelectrons
#			    PE = numpy.random.normal(anode_current*T/e, noise_current*T/e, None)
#
#			    # A/D converter : Current --> ADC counts
#			    ADC = self.get_ADC_value(pixel, PE)

			if (self.configs.detector_mode == 'Pulse') :
			    # expectation
			    E = round(QE*photon_flux*T)

			    # get photoelectron-signal
			    signal = numpy.random.poisson(E, None)

			    # get dark count
			    dark  = 0
			    #dark  = numpy.random.normal(0, round(D*T), None)

			    # get readout noise
			    noise = 0
			    #noise = numpy.random.normal(0, NA, None)

			    # get total pulses
			    pulse = (signal + background) + dark + noise

			    # A/D converter : Pulse --> ADC counts
			    ADC = self.get_ADC_value(pixel, pulse)

			if (self.configs.spatiocyte_bc_switch == True) :
			    cell_pixel[i][j] = ADC*bc[i][j]
			else :
			    cell_pixel[i][j] = ADC


                # Background (No photon signal)
                camera_pixel = numpy.zeros(shape=(Nw_pixel, Nh_pixel))

                for i in range(Nw_pixel) :
                    for j in range(Nh_pixel) :

			# set pixel position
			pixel = (i, j)

			photon_flux, background = 0, 0

#			if (self.configs.detector_mode == 'Current') :
#
#			    # get current at cathode (photoelectric current)
#			    e = self.configs.electron_charge
#			    cathode_current = e*(signal + background)
#
#			    # get noise
#			    noise_current = self.get_noise_analog(cathode_current + background)
#
#			    # get anode current (Gain x cathode current)
#			    G = self.configs.detector_emgain
#			    anode_current = G*(cathode_current + background)
#
#			    # get total charge
#			    PE = numpy.random.normal(anode_current*T/e, noise_current*T/e, None)
#
#			    # A/D converter : Charge --> ADC count
#			    ADC = self.get_ADC_value(pixel, PE)

			if (self.configs.detector_mode == 'Pulse') :
			    # expectation
			    E = round(QE*photon_flux*T)

			    # get photoelectron-signal
			    signal = numpy.random.poisson(E, None)

			    # get dark count
			    dark = 0
			    #dark = numpy.random.normal(0, round(D*T), None)

			    # get noise
			    noise = 0
			    #noise = numpy.random.normal(0, NA, None)

			    # get pulse
			    pulse = (signal + background) + dark + noise

			    # A/D converter : Pulse  --> ADC count
			    ADC = self.get_ADC_value(pixel, pulse)

			camera_pixel[i][j] = ADC


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


def bounded(value, lower_bound=None, upper_bound=None):
    if lower_bound is not None and value < lower_bound:
        return lower_bound
    if upper_bound is not None and value > upper_bound:
        return upper_bound
    return value
