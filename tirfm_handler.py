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

import tirfm_settings

from scipy.special import j0
from PIL import Image

IMAGE_SIZE_LIMIT=3000


class VisualizerError(Exception):

    "Exception class for visualizer"

    def __init__(self, info):
        self.__info = info

    def __repr__(self):
        return self.__info

    def __str__(self):
        return self.__info



#class TIRFMSettings(Settings) :
class TIRFMSettings() :

    '''
    TIRFM Visualization setting class

	Incident Beam
	Fluorophore
	Mirror
	Dichroic Mirror
	Excitation Filter
	Emission Filter
	Pinhole
	Tube Lenz
	Camera
    '''

    def __init__(self, user_settings_dict = None):

        # default setting
        settings_dict = tirfm_settings.__dict__.copy()
        #settings_dict_tirfm = tirfm_settings.__dict__.copy()
        #settings_dict.update(settings_dict_tirfm)

        # user setting
        if user_settings_dict is not None:
            if type(user_settings_dict) != type({}):
                print 'Illegal argument type for constructor of Settings class'
                sys.exit()
            settings_dict.update(user_settings_dict)

        for key, val in settings_dict.items():
            if key[0] != '_': # Data skip for private variables in setting_dict.
                if type(val) == type({}) or type(val) == type([]):
                    copy_val = copy.deepcopy(val)
                else:
                    copy_val = val
                setattr(self, key, copy_val)


    def _set_data(self, key, val) :
    
        if val != None:
            setattr(self, key, val)



    def set_Particles(self, p_list = None) :

	print '--- Particles List :'

        self._set_data('particles_list', p_list)

	print p_list



    def set_IncidentBeam(self,  wave_length = None,
                                intensity = None,
                                pinhole_radius = None,
                                excitation = None ) :

        print '--- Incident Beam Condition :'

        self._set_data('beam_switch', True)
        self._set_data('beam_wavelength', wave_length)
        self._set_data('beam_intensity', intensity)
        self._set_data('beam_pinhole_radius', pinhole_radius)

        print '\tWave Length = ', self.beam_wavelength, 'nm'
        print '\tIntensity   = ', self.beam_intensity, 'W/cm^2 (=joule/sec cm^2)'
        print '\tPinhole Radius = ', self.beam_pinhole_radius, 'nm'


        if (excitation == None) :
            print '\tExcitation Filter OFF'
        else :
            print '\tExcitation Filter ON'

            filename = './catalog/excitation/' + excitation + '.csv'

            try:
                csvfile = open(filename)
                lines = csvfile.readlines()

                header = lines[0:5]
                data   = lines[6:]

                excitation_header = []
                excitation_filter = []

                for i in range(len(header)) :
                    dummy  = header[i].split('\r\n')
                    a_data = dummy[0].split(',')
                    excitation_header.append(a_data)
                    print '\t', a_data

                for i in range(len(data)) :
                    dummy0 = data[i].split('\r\n')
                    a_data = dummy0[0].split(',')
                    excitation_filter.append(a_data)

            except Exception:
                print 'Error : ', filename, ' is NOT found'
                exit()

            ####
            self.excitation_eff = self.set_efficiency(excitation_filter)
            self._set_data('excitation_switch', True)



    def set_Fluorophore(self, fluorophore_type = None,
				wave_length = None,
				width = None,
				cutoff = None,
                        	file_name_format = None ) :

	if (fluorophore_type == 'Gaussian') :

            print '--- Fluorophore : Point Spreading Function [%s]' % (fluorophore_type)

            self._set_data('fluorophore_type', fluorophore_type)
            self._set_data('psf_wavelength', wave_length)
            self._set_data('psf_width', width)
            self._set_data('psf_cutoff', cutoff)
            self._set_data('psf_file_name_format', file_name_format)

	    index = self.psf_wavelength - self.fluorophore_wavelength[0]
	    self.fluoex_eff[index] = 100
	    self.fluoem_eff[index] = 100

            print '\tWave Length   = ', self.psf_wavelength, 'nm'
            print '\tRadial  Width = ', self.psf_width[0], 'nm'
            print '\tLateral Width = ', self.psf_width[1], 'nm'
            print '\tRadial  Cutoff = ', self.psf_cutoff[0], 'nm'
            print '\tLateral Cutoff = ', self.psf_cutoff[1], 'nm'



	elif (fluorophore_type == 'Point-like') :

            print '--- Fluorophore : Point Spreading Function [%s]' % (fluorophore_type)

            self._set_data('fluorophore_type', fluorophore_type)
            self._set_data('psf_wavelength', wave_length)
            self._set_data('psf_file_name_format', file_name_format)

            index = self.psf_wavelength - self.fluorophore_wavelength[0]
            self.fluoex_eff[index] = 100
            self.fluoem_eff[index] = 100

            print '\tWave Length   = ', self.psf_wavelength, 'nm'



	else :

	    print '--- Fluorophore :'

            filename = './catalog/fluorophore/' + fluorophore_type + '.csv'

            try:
		csvfile = open(filename)
		lines = csvfile.readlines()

		header = lines[0:5]
		data   = lines[5:]

		fluorophore_header = []
		fluorophore_excitation = []
		fluorophore_emission = []

		for i in range(len(header)) :
                    dummy  = header[i].split('\r\n')
                    a_data = dummy[0].split(',')
                    fluorophore_header.append(a_data)
		    print '\t', a_data

		for i in range(len(data)) :
                    dummy0 = data[i].split('\r\n')
                    a_data = dummy0[0].split(',')
	
		    if   (len(a_data) == 1 and a_data[0] == 'Excitation') : flag = 0
		    elif (len(a_data) == 1 and a_data[0] == 'Emission'  ) : flag = 1
		    else :
		        if (flag == 0) : 
			    fluorophore_excitation.append(a_data)
		        else : 
			    fluorophore_emission.append(a_data)


            except Exception:
                print 'Error : ', filename, ' is NOT found'
                exit()

	    ####
            self._set_data('fluorophore_type', fluorophore_type)
	    self.fluoex_eff = self.set_efficiency(fluorophore_excitation)
	    self.fluoem_eff = self.set_efficiency(fluorophore_emission)

	# Normalization
	norm = sum(self.fluoex_eff)
	self.fluoex_norm = self.fluoex_eff/norm

	norm = sum(self.fluoem_eff)
	self.fluoem_norm = self.fluoem_eff/norm



    def set_Mirror(self, position = None ) :

        print '--- Mirror :'

        self._set_data('mirror_position', position)

        print '\tMirror Position = ', self.mirror_position



    def set_Objective(self, 
			NA = None,
			Ng = None,
			Nm = None,
			mag = None,
			efficiency = None,
			thickness  = None
			) :

	print '--- Objective :'

        self._set_data('objective_switch', True)
        self._set_data('objective_NA', NA)
        self._set_data('objective_Ng', Ng)
        self._set_data('objective_Nm', Nm)
        self._set_data('objective_mag', mag)
        self._set_data('objective_efficiency', efficiency)
        self._set_data('objective_glassthickness', thickness)
	self._set_angle()

	print '\tNA = ', self.objective_NA
	print '\tN(glass)  = ', self.objective_Ng
	print '\tN(medium) = ', self.objective_Nm
	print '\tMagnification = ', self.objective_mag, 'x'
	print '\tTransmission Efficiency = ', self.objective_efficiency
	print '\tsin(max angle)      = ', self.objective_sin_max
	print '\tsin(critical angle) = ', self.objective_sin_critical



    def _set_angle(self) :

	self.objective_sin_max = self.objective_NA/self.objective_Ng
	self.objective_sin_critical = self.objective_Nm/self.objective_Ng



    def set_DichroicMirror(self, dm = None) :

        print '--- Dichroic Mirror :'

	filename = './catalog/dichroic/' + dm + '.csv'

	try:
	    csvfile = open(filename)
	    lines = csvfile.readlines()

	    header = lines[0:5]
	    data   = lines[6:]

	    dichroic_header = []
	    dichroic_mirror = []

	    for i in range(len(header)) :
		dummy  = header[i].split('\r\n')
		a_data = dummy[0].split(',')
		dichroic_header.append(a_data)
		print '\t', a_data

	    for i in range(len(data)) :
		dummy0 = data[i].split('\r\n')
		a_data = dummy0[0].split(',')
		dichroic_mirror.append(a_data)

        except Exception:
            print 'Error : ', filename, ' is NOT found'
	    exit()


        self.dichroic_eff = self.set_efficiency(dichroic_mirror)

	self._set_data('dichroic_switch', True)



    def set_EmissionFilter(self, emission = None) :

        print '--- Emission Filter :'

        filename = './catalog/emission/' + emission + '.csv'

        try:
            csvfile = open(filename)
            lines = csvfile.readlines()

            header = lines[0:5]
            data   = lines[6:]

	    emission_header = []
	    emission_filter = []

            for i in range(len(header)) :
                dummy  = header[i].split('\r\n')
                a_data = dummy[0].split(',')
                emission_header.append(a_data)
		print '\t', a_data

            for i in range(len(data)) :
                dummy0 = data[i].split('\r\n')
                a_data = dummy0[0].split(',')
                emission_filter.append(a_data)

        except Exception:
            print 'Error : ', filename, ' is NOT found'
            exit()

        self.emission_eff = self.set_efficiency(emission_filter)
        self._set_data('emission_switch', True)



    def set_TubeLens(self, mag = None, wl_range = None, efficiency = None) :

        print '--- Tube Lens :'
	
        self._set_data('tubelens_switch', True)
        self._set_data('tubelens_mag', mag)
        self._set_data('tubelens_range', wl_range)
        self._set_data('tubelens_efficiency', efficiency)

        print '\tMagnification = ', self.tubelens_mag, 'x'
        print '\tWavelength Range = ', self.tubelens_range, 'nm'
        print '\tTransmission Efficiency = ', self.tubelens_efficiency



    def set_Pinhole(self, radius) :

	print '--- Pinhole :'

        self._set_data('pinhole_radius', radius)

        print '\tRadius = ', self.pinhole_radius, 'nm'


    def set_Camera(self, camera = None,
		   image_size = None,
                   pixel_length= None,
                   focal_point = None,
                   base_position = None,
                   zoom = None,
                   start_time = None,
                   end_time = None,
                   fps = None,
                   exposure_time = None,
		   sat_charge = None,
		   ADC_bit = None,
		   ADC_const = None,
		   ADC_offset = None,
		   dark_current = None,
		   readout = None,
		   excess = None,
		   emgain = None
                   ):

        print '--- Camera :'

        
        try:
	    filename = './catalog/camera/' + camera + '.csv'

            csvfile = open(filename)
            lines = csvfile.readlines()

            header = lines[0:12]
            data   = lines[13:]

            camera_header = []
            camera_QEdata = []

            for i in range(len(header)) :
                dummy  = header[i].split('\r\n')
                a_data = dummy[0].split(',')
                camera_header.append(a_data)
                #print '\t', a_data

	    image_size	 = (int(camera_header[5][1]), int(camera_header[5][1]))
	    pixel_length = float(camera_header[6][1])
	    sat_charge   = float(camera_header[7][1])
	    readout	 = float(camera_header[8][1])
	    dark_current = float(camera_header[9][1])
	    excess	 = float(camera_header[10][1])
	    emgain	 = float(camera_header[11][1])

            for i in range(len(data)) :
                dummy0 = data[i].split('\r\n')
                a_data = dummy0[0].split(',')
                camera_QEdata.append(a_data)

	    self.camera_blue  = self.set_efficiency(camera_QEdata, 1)
	    self.camera_green = self.set_efficiency(camera_QEdata, 2)
	    self.camera_red   = self.set_efficiency(camera_QEdata, 3)



        except Exception:

	    print 'Error : ', filename, ' file is NOT found'
	    print '\tUse Perfect camera (Noise free)\n'

	    camera = 'Perfect'
            #exit()

        self._set_data('camera_switch', True)
        self._set_data('camera_type', camera)
        self._set_data('camera_image_size', image_size)
        self._set_data('camera_pixel_length', pixel_length)
        self._set_data('camera_focal_point', focal_point)
        self._set_data('camera_base_position', base_position)
        self._set_data('camera_zoom', zoom)
        self._set_data('camera_start_time', start_time)
        self._set_data('camera_end_time', end_time)
        self._set_data('camera_fps', fps)
        self._set_data('camera_exposure_time', exposure_time)
        self._set_data('camera_sat_charge', sat_charge)
        self._set_data('camera_ADC_bit', ADC_bit)
        self._set_data('camera_ADC_const', ADC_const)
        self._set_data('camera_ADC_offset', ADC_offset)
        self._set_data('camera_dark_current', dark_current)
        self._set_data('camera_readout', readout)
        self._set_data('camera_excess', excess)
        self._set_data('camera_emgain', emgain)

	print '\tCamera Type : ', self.camera_type
        print '\tImage Size  = ', self.camera_image_size[0], 'x', self.camera_image_size[1]
        print '\tPixel Size  = ', self.camera_pixel_length, 'm/pixel'
        print '\tFocal Point = ', self.camera_focal_point
        print '\tPosition    = ', self.camera_base_position
        print '\tZoom	     = ', self.camera_zoom
        print '\tStart Time  = ', self.camera_start_time, 'sec'
        print '\tEnd   Time  = ', self.camera_end_time, 'sec'
        print '\tFrame Rate  = ', self.camera_fps, 'frames/sec'
        print '\tExposure Time = ', self.camera_exposure_time, 'sec'
        print '\tSat. Charge   = ', self.camera_sat_charge, 'electron'
        print '\tA/D Converter = ', self.camera_ADC_bit, 'bit'
        print '\tADC Const  = ', self.camera_ADC_const, 'ADC count/electron'
        print '\tADC Offset = ', self.camera_ADC_offset, 'ADC count'
        print '\tDark Current = ', self.camera_dark_current
        print '\tReadout = ', self.camera_readout
        print '\tExcess	 = ', self.camera_excess
        print '\tEM gain = ', self.camera_emgain



    def set_Movie(self, image_file_dir = None,
			movie_filename='./movies/movie.mp4', \
			cleanup_image_file_dir=False) :

	if image_file_dir is None:
	    image_file_dir = tempfile.mkdtemp(dir=os.getcwd())
	    cleanup_image_file_dir = True
            
	self._set_data('movie_image_file_dir', image_file_dir)
	self._set_data('movie_cleanup_image_file_dir', cleanup_image_file_dir)
	self._set_data('movie_filename', movie_filename)



    def set_DataFile(self, hdf5_file_path_list) :

	# read hdf5 lattice file
	for hdf5_file_path in hdf5_file_path_list :

	    try :

		hdf5_file = h5py.File(hdf5_file_path, 'r')
	
		species = hdf5_file['species']
		lattice = hdf5_file['lattice_info/HCP_group']
	        dataset = hdf5_file['data']
	
		### particle data in time-series
		data = []

		for i in dataset :

		    data_i= dataset[i]
	            time = data_i.attrs['t']
		    particles = hdf5_file['data/'+str(i)+'/particles']
		    element = [time, particles]
		    data.append(element)

		data.sort(lambda x, y:cmp(x[0], y[0]))

		# get data
		self._set_data('spatiocyte_data', data)

		# get species properties
		self._set_data('spatiocyte_species_id', map(lambda x : x[0], species))
		self._set_data('spatiocyte_index',      map(lambda x : x[1], species))
		self._set_data('spatiocyte_diffusion',  map(lambda x : x[3], species))
		self._set_data('spatiocyte_radius',     map(lambda x : x[2], species))

		# get lattice properties
		self._set_data('spatiocyte_lattice_id', map(lambda x : x[0], lattice))
		self._set_data('spatiocyte_lengths',    map(lambda x : x[1], lattice))
		self._set_data('spatiocyte_VoxelRadius',   lattice[0][2])
		self._set_data('spatiocyte_theNormalizedVoxelRadius', lattice[0][3])
		self._set_data('spatiocyte_theStartCoord', lattice[0][4])
		self._set_data('spatiocyte_theRowSize',    lattice[0][6])
		self._set_data('spatiocyte_theLayerSize',  lattice[0][5])
		self._set_data('spatiocyte_theColSize',    lattice[0][7])

		#hdf5_file.close()

            except Exception, e :
	                if not self.ignore_open_errors:
	                    raise
	                print 'Ignoring error: ', e


	# Visualization error	
	if self.spatiocyte_species_id is None:
	    raise VisualizerError('Cannot find species_id in any given hdf5 files')

	if len(self.spatiocyte_data) == 0:
	    raise VisualizerError('Cannot find particles dataset in any given hdf5 files: ' \
	                    + ', '.join(hdf5_file_path_list))

	if len(self.spatiocyte_index) == 0 :
	    raise VisualizerError('Cannot find lattice dataset in any given hdf5 files: ' \
	                    + ', '.join(hdf5_file_path_list))



    def set_efficiency(self, array, index=1) :

	if (len(array[0]) < 3) : index = 1

        N = len(self.fluorophore_wavelength)
        efficiency = numpy.array([0.0 for i in range(N)])

        for i in range(N) :
            wl = self.fluorophore_wavelength[i]

            for j in range(len(array)) :

                length = float(array[j][0])
                eff = float(array[j][index])

                if (length/wl == 1) :
                    efficiency[i] = eff


        return efficiency



    def set_depth(self, Magnitude) :

	wave_length = self.beam_wavelength

        n1 = self.objective_Ng
        n2 = self.objective_Nm

        length = math.sqrt(1 - self.objective_sin_max**2)/self.objective_sin_max
        tan_th = self.mirror_position/length
        sin2   = tan_th/math.sqrt(1 + tan_th**2)

        ref = n1**2*sin2 - n2**2

        if (ref > 0) :

            self.penetration_depth = wave_length/(4*math.pi*math.sqrt(ref))

            print '--- TIRF Configuration : '

        else :

            self.penetration_depth = float('inf')

            print '--- EPIF Configuration : '

	# focal depth
	depth_1st = self.objective_Ng/(2.0*self.objective_NA**2)*wave_length
	depth_2nd = self.objective_Ng/(7*Magnitude*self.objective_NA)

	self.focal_depth = depth_1st + depth_2nd


	print '\tPenetration Depth = ', self.penetration_depth, 'nm'
	print '\tFocal Depth = ', self.focal_depth, 'nm'



    def get_Photons(self, z, wave_length) :

        # plank const * speed of light
        hc = 2.00e-25 # J m

        # single photon energy
        E_wl = hc/(wave_length*1.0e-9)
        area = numpy.pi*self.spatiocyte_VoxelRadius**2

        # calculate the number of photons per voxel
        N0 = (self.beam_intensity*1.0e+4/E_wl)*area*self.camera_exposure_time

        N_pd = numpy.exp(-z/self.penetration_depth)
        N_fd = numpy.exp(-0.5*(z/self.focal_depth)**2)

        #N_photons_in = N0*N_pd*N_fd
        N_photons_in =  numpy.array(map(lambda x : N0*x, N_pd*N_fd))

        #################################################################
        #
        # Note : Adding Quantum mechanics of Photons to Molecules
        #
        #       Using random number generator to simulate
        #       the energy transition of molecular vibrational states
        #
        #################################################################

        N_photons_out = 1.0e-6*N_photons_in

        return N_photons_out



    def set_PSF(self) :

        wave_length = self.fluorophore_wavelength

	r = self.fluorophore_radial # in nm-scale
	z = self.fluorophore_depth  # in nm-scale

        # Fluorophores Emission Intensity (wave_length)
        I = self.fluoem_norm

        # Photon Transmission Efficiency
        if (self.dichroic_switch == True) :
            I = I*0.01*self.dichroic_eff

        if (self.emission_switch == True) :
            I = I*0.01*self.emission_eff


	# count photons
	N_ph = map(lambda x : I*x, self.get_Photons(z, wave_length))

	# Detector : Quantum Efficiency
	self.fluorophore_rgb[:,2] = map(lambda x : sum(x), map(lambda x : x*self.camera_blue,  N_ph))
	self.fluorophore_rgb[:,1] = map(lambda x : sum(x), map(lambda x : x*self.camera_green, N_ph))
	self.fluorophore_rgb[:,0] = map(lambda x : sum(x), map(lambda x : x*self.camera_red,   N_ph))

	# count photoelectrons
	N_b = sum(map(lambda x : x*self.camera_blue,  N_ph))
	N_g = sum(map(lambda x : x*self.camera_green, N_ph))
	N_r = sum(map(lambda x : x*self.camera_red,   N_ph))

	# for Normalization
	norm = map(lambda x : True if x > 1e-2 else False, N_ph[0])

	# Point Spread Function
	if (self.fluorophore_type == 'Gaussian') :

            I0 = 1.0
            Ir = sum(map(lambda x : x*numpy.exp(-0.5*(r/self.psf_width[0])**2), Norm_array))
            Iz = sum(map(lambda x : x*numpy.exp(-0.5*(z/self.psf_width[1])**2), Norm_array))

	    self.fluorophore_psf = numpy.array(map(lambda x : I0*Ir*x, Iz))


        elif (self.fluorophore_type == 'Point-like') :

            I0 = 1.0
            Ir = sum(map(lambda x : x*numpy.array(map(lambda x : 1.00 if x == 0 else 0.00, r)), Norm_array))
            Iz = sum(map(lambda x : x*numpy.array(map(lambda x : 1.00 if x == 0 else 0.00, z)), Norm_array))

            self.fluorophore_psf = numpy.array(map(lambda x : I0*Ir*x, Iz))


	else :

	    # make the norm and wave_length array shorter
	    re_wavelength = []

	    for i in range(len(norm)) :

		if norm[i] is True :
		    re_wavelength.append(wave_length[i])

	    self.get_PSF(r, z, numpy.array(re_wavelength))

#	    for k in range(len(self.fluorophore_wavelength)) :
#
#		wave_length = self.fluorophore_wavelength[k]
#
#		if (Norm_array[k] is False) :
#		     continue
#
#		print wave_length, 'nm'
#
#		self.get_PSF(r, z, wave_length)



	# Normalization
	self.fluorophore_psf = self.fluorophore_psf/sum(norm)



    def get_PSF(self, r, z, wave_length) :

	NA = self.objective_NA
	N = 50
	drho = 1.0/N
	rho = numpy.array([(i+1)*drho for i in range(N)])

	k = 2.0*numpy.pi/wave_length
	alpha = k*NA
	gamma = k*(NA/2)**2

	J0 = numpy.array(map(lambda y : map(lambda x : j0(x*y*rho), r), alpha))
	Y  = numpy.array(map(lambda y : map(lambda x : 2*numpy.exp(-2*1.j*x*y*rho**2)*rho*drho, z), gamma))
	#I = numpy.array(map(lambda c0, c1 : map(lambda b : map(lambda a : 2*j0(a*c0*rho)*numpy.exp(-2*1.j*b*c1*rho**2)*rho*drho, r), z), alpha, gamma))
	I  = numpy.array([numpy.array(map(lambda x : x*J0[i], Y[i])) for i in range(len(wave_length))])
	#I.sum(axis=3)
	I_abs = map(lambda x : abs(x)**2, I.sum(axis=3))
	self.fluorophore_psf = sum(I_abs)

	print self.fluorophore_psf[0]

#	for i in range(len(wave_length)) :
#
#	    I  = numpy.array(map(lambda x : x*J0[i], Y[i]))
#	    I_sum = I.sum(axis=2)
#	    I_abs = map(lambda x : abs(x)**2, I_sum)
#
#	    self.fluorophore_psf += I_abs
#
#	    print wave_length[i], I_abs[0][0]

#       J0 = map(lambda x : j0(x*alpha*rho), r)
#       Y  = map(lambda x : 2*numpy.exp(-2*1.j*x*gamma*rho**2)*rho*drho, z)
#       I  = numpy.array(map(lambda x : x*J0, Y))
#
#	I_sum = I.sum(axis=2)
#	I_abs = map(lambda x : abs(x)**2, I_sum)
#
#	self.fluorophore_psf += I_abs



#    def get_PSF(self, r, z, wave_length) :
#
#	NA = self.objective_NA
#        N = 100
#        drho = 1.0/N
#        rho = numpy.array([i*drho for i in range(N)])
#
#        k = 2.0*numpy.pi/wave_length
#	alpha = k*NA
#	gamma = k*(NA/2)**2
#
#	J0 = map(lambda x : j0(x*alpha*rho), r)
#	Y  = map(lambda x : 2*numpy.exp(-2*1.j*x*gamma*rho**2)*rho*drho, z)
#
#	I  = numpy.array(map(lambda x : x*J0, Y))
#
#	I_sum = I.sum(axis=2)
#	I_abs = map(lambda x : abs(x)**2, I_sum)
#
#	self.fluorophore_psf += I_abs
#
#	print I_abs[0][0]



    def get_Noise(self, signal) :

	Nr = self.camera_readout
	DC = self.camera_dark_current
	Fn = numpy.sqrt(self.camera_excess)
	M  = self.camera_emgain

	# get noise -- defined by HAMAMATSU Photonics
	sigma2 = Fn**2*(signal + DC*self.camera_exposure_time)+ (Nr/M)**2
	noise  = numpy.sqrt(sigma2)

	return noise



    def set_SNR(self) :

	for j in range(len(self.fluorophore_wavelength)) :

	    index = self.fluorophore_wavelength[j] - self.fluorophore_wavelength[0]

            QEff = self.camera_green[index]

	    for i in range(len(self.photon_number)) :
		# signal
		signal = self.photon_number[i]

		# get SNR and relative SNR
		self.absolute_snr[j][i] = (QEff*signal)/self.get_Noise(QEff*signal)
		self.relative_snr[j][i] = self.absolute_snr[j][i]/self.ideal_snr[i]



#class TIRFMVisualizer(Visualizer) :
class TIRFMVisualizer() :

	'''
	TIRFM Visualization class of e-cell simulator
	'''

	def __init__(self, settings=TIRFMSettings()) :

		assert isinstance(settings, TIRFMSettings)
		self.settings = settings

		"""
		Check and create the folder for image file.
		"""
		if not os.path.exists(self.settings.movie_image_file_dir):
		    os.makedirs(self.settings.movie_image_file_dir)
		else:
		    for file in os.listdir(self.settings.movie_image_file_dir):
			os.remove(os.path.join(self.settings.movie_image_file_dir, file))

                """
                Scaling Factor
                """
                voxel_radius = self.settings.spatiocyte_VoxelRadius

                view = self.settings.camera_pixel_length/(2.0*voxel_radius)
                Mag  = self.settings.objective_mag*self.settings.tubelens_mag
                zoom = self.settings.camera_zoom

                self.image_scaling = (view/Mag)*zoom

                """
                Focal and Penetration Depth
                """
                self.settings.set_depth(Mag)

                """
                Point Spread Function (PSF)
                """
                self.settings.set_PSF()

                """
                Image Size and Boundary
                """
                width  = int(self.settings.camera_image_size[0])
                height = int(self.settings.camera_image_size[1])

                if width > IMAGE_SIZE_LIMIT or height > IMAGE_SIZE_LIMIT :
                        raise VisualizerErrror('Image size is bigger than the limit size')

                # camera's focal position
		focal = numpy.array(self.settings.camera_focal_point)

		# get length
                pixel_z = self.settings.spatiocyte_lengths[0][2]/self.image_scaling
                pixel_y = self.settings.spatiocyte_lengths[0][1]/self.image_scaling
                pixel_x = self.settings.spatiocyte_lengths[0][0]/self.image_scaling

		x_min = 0
		x_max = int(pixel_x)

                y_min = int(width/2.0  - pixel_y*focal[1] - 1)
                y_max = int(height/2.0 + pixel_y*focal[1] + 1)

                z_min = int(width/2.0  - pixel_z*focal[2] - 1)
                z_max = int(height/2.0 + pixel_z*focal[2] + 1)

                if (z_min < 0) : z_min = 0
                if (z_max > width)  : z_max = width
                #if (x_min < 0) : x_min = 0
                #if (x_max > width)  : x_max = width

                if (y_min < 0) : y_min = 0
                if (y_max > height) : y_max = height

		# image size
		self.img_width  = width
		self.img_height = height

                # image max/min position
                self.img_min = (x_min, y_min, z_min)
                self.img_max = (x_max, y_max, z_max)

		print self.img_min
		print self.img_max



        def __del__(self):

            if self.settings.movie_cleanup_image_file_dir :

                for parent_dir, dirs, files in os.walk(self.settings.movie_image_file_dir, False) :
                    for file in files :
                        os.remove(os.path.join(parent_dir, file))

                    os.rmdir(parent_dir)



	def _get_coordinate(self, aCoord) :

	        """
		get (column, layer, row) coordinate
	        """
		start_coord = self.settings.spatiocyte_theStartCoord
		row_size    = self.settings.spatiocyte_theRowSize
		layer_size  = self.settings.spatiocyte_theLayerSize
		col_size    = self.settings.spatiocyte_theColSize

                aGlobalCol   =  (aCoord-start_coord)/(row_size*layer_size)
                aGlobalLayer = ((aCoord-start_coord)%(row_size*layer_size))/row_size
                aGlobalRow   = ((aCoord-start_coord)%(row_size*layer_size))%row_size

                """
		get (x, y, z) coordinate
                """
		norm_voxel_radius = self.settings.spatiocyte_theNormalizedVoxelRadius

                theHCPk = norm_voxel_radius/math.sqrt(3.0)
                theHCPh = norm_voxel_radius*math.sqrt(8.0/3.0)
                theHCPl = norm_voxel_radius*math.sqrt(3.0)

	        point_y = (aGlobalCol%2)*theHCPk + theHCPl*aGlobalLayer
	        point_z = aGlobalRow*2*norm_voxel_radius + ((aGlobalLayer+aGlobalCol)%2)*norm_voxel_radius
	        point_x = aGlobalCol*theHCPh

		return point_x, point_y, point_z



	def _get_center(self) :

		norm_voxel_radius = self.settings.spatiocyte_theNormalizedVoxelRadius

                theHCPk = norm_voxel_radius/math.sqrt(3.0)
                theHCPh = norm_voxel_radius*math.sqrt(8.0/3.0)
                theHCPl = norm_voxel_radius*math.sqrt(3.0)

		lengths_z = self.settings.spatiocyte_lengths[0][2]
		lengths_y = self.settings.spatiocyte_lengths[0][1]
		lengths_x = self.settings.spatiocyte_lengths[0][0]

		row0   = lengths_z/2.0 + 4.0*norm_voxel_radius
		layer0 = lengths_y/2.0 + 2.0*theHCPl
		col0   = lengths_x/2.0 + 2.0*theHCPh

		return row0, layer0, col0


	def get_position(self, pos) :

                # normal vector
                norm  = map(operator.sub, self.settings.camera_base_position, self.settings.camera_focal_point)
                norm_len2 = norm[0]**2 + norm[1]**2 + norm[2]**2
                norm_len = math.sqrt(norm_len2)

                # camera's focal position
                focal = numpy.array(self.settings.camera_focal_point)

                coef_d = -norm[0]*focal[0]-norm[1]*focal[1]-norm[2]*focal[2]
		par_xy = [0, self.img_width/2, self.img_height/2]

                # transform matrix
                tmat = trans_mat(norm)

		# rotate
		p_inp = norm[0]*pos[0]+norm[1]*pos[1]+norm[2]*pos[2]+coef_d
		dis = math.fabs(p_inp)/norm_len

		# calculate foot of perpendicular
		t0 = -p_inp/norm_len2
		pos_ft =[norm[0]*t0+pos[0], norm[1]*t0+pos[1], norm[2]*t0+pos[2]]

		# transform the foot coordinate to x-y plane
		pos_xy = numpy.dot(tmat, numpy.array(pos_ft - focal))
		pos_xyp = pos_xy + par_xy

		pos = (pos_xyp[0,0], pos_xyp[0,1], pos_xyp[0,2])

		return pos



	def _Data2Frame(self) :

        	frame_data = []
        	frame_time = []
        	expos_time = []

		frame_interval = 1.0/self.settings.camera_fps

        	# set the frame and exposure time (start/end time)
        	counter = 1
        	ignore_dtime = frame_interval/1.0e+5

        	while True :

        	    ft = self.settings.camera_start_time + frame_interval*counter
        	    frame_time.append(ft)
        	    et = ft - self.settings.camera_exposure_time
        	    if (et < 0.0) : et = 0.0
        	    expos_time.append(et)
        	    if (ft >= self.settings.camera_end_time) : break

        	    counter += 1


        	# modify step-to-step data format to frame-to-frame data format
        	for step in range(len(frame_time)) :

        	    ft = frame_time[step]
        	    et = expos_time[step]

        	    frame_elem = [ft]
        	    last_index = 0

        	    for index in range(len(self.settings.spatiocyte_data)):

        	        if index == 0 : continue

        	        st = self.settings.spatiocyte_data[index][0]

        	        if(et - st <= ignore_dtime and st - ft <= ignore_dtime) :

        	            st_f = self.settings.spatiocyte_data[index-1][0]
        	            stay_time  = min(st - st_f, st - et)
        	            norm_stime = stay_time/self.settings.camera_exposure_time
	
			    element = (norm_stime, self.settings.spatiocyte_data[index][1])
        	            frame_elem.append(element)

        	            last_index = index

        	    # check last data
        	    if last_index == 0 : continue

        	    st = self.settings.spatiocyte_data[last_index][0]
        	    stay_time = ft - st

        	    if stay_time > ignore_dtime :

        	        norm_stime = stay_time/self.settings.camera_exposure_time
        	        element = (norm_stime, self.settings.spatiocyte_data[last_index][1])
        	        frame_elem.append(element)
	
        	    frame_data.append(frame_elem)


        	return frame_data



	def _get_Signal(self, p0, pixel) :

		# get signal (photoelectron counting)
		signal = numpy.array([0, 0, 0])

		# EM gain
		gain = self.settings.camera_emgain

		x, z, y = pixel

		for i in range(len(p0)) :

		    x_i = p0[i][0]
		    y_i = p0[i][1]
		    z_i = p0[i][2]

		    r = numpy.sqrt((z - z_i)**2 + (y - y_i)**2)
		    d = abs(x - x_i)

		    # convert pixel to real(nm) scale
		    voxel_radius = self.settings.spatiocyte_VoxelRadius

                    r = r*(2.0*voxel_radius*self.image_scaling/1e-9)
                    d = d*(2.0*voxel_radius*self.image_scaling/1e-9)

		    # get signal
		    if (int(d) < len(self.settings.fluorophore_depth) and
		        int(r) < len(self.settings.fluorophore_radial) and
		        int(r) < self.settings.pinhole_radius) :

			N_pe = self.settings.fluorophore_rgb[int(d)]

		        signal[0] += N_pe[0]*gain*self.settings.fluorophore_psf[int(d)][int(r)]
		        signal[1] += N_pe[1]*gain*self.settings.fluorophore_psf[int(d)][int(r)]
		        signal[2] += N_pe[2]*gain*self.settings.fluorophore_psf[int(d)][int(r)]

		# get noise
		noise = [self.settings.get_Noise(signal[i]) for i in range(3)]

		# random generator
		#I_signal = numpy.random.poisson(signal, None)
		#I_noise  = numpy.random.poisson(noise,  None)

		#I = I_signal + I_noise
		#I = signal

                # convert photoelectron to ADC values(Grayscale)
                k = self.settings.camera_ADC_const
                ADC0 = self.settings.camera_ADC_offset
                ADC  = k*signal + ADC0

                # Rescale to 8-bit
                ADC = ADC*(2.0**8-1)/(2.0**self.settings.camera_ADC_bit-1)

		return (int(ADC[0]), int(ADC[1]), int(ADC[2]))



	def output_frames(self, num_div=1) :
		"""
	        Output Images
	        """
                # - Frame2Frame data set
                # convert step-to-step dataset to frame-to-frame one
                frame_data = self._Data2Frame()

	        # create frame data composed by frame element data
	        for i in range(len(frame_data)) :

		    # set image file name
		    image_file_name = os.path.join(self.settings.movie_image_file_dir, \
						self.settings.movie_image_file_name_format % i)

		    # initialize tirfm image
		    tirfm_image = Image.new("RGB", (self.img_width, self.img_height), (0,0,0))

		    #####
                    time  = frame_data[i][0]

                    print 'time : ', time, ' sec'

		    p0 = []

                    for j in range(1, len(frame_data[i])-1) :

			n_st = frame_data[i][j][0]
                        c_id = map(lambda x : x[0], frame_data[i][j][1])
                        s_id = map(lambda x : x[1], frame_data[i][j][1])
                        l_id = map(lambda x : x[2], frame_data[i][j][1])

			for k in range(len(c_id)) :
			  # Raf(Membrane)
			  #if (s_id[k] != 34) :
			  # ERK-pp
			  #if (s_id[k] == 3 or s_id[k] == 41) :
			  #if (s_id[k] == 6 or s_id[k] == 12) :
			  # ERK-all
			  #if (s_id[k] != 37 or s_id[k] != 38 or s_id[k] != 44) :
			  #if (s_id[k] != 0 or s_id[k] != 1 or s_id[k] != 15) :

			    # particles coordinate in real(nm) scale
                            x, y, z = self._get_coordinate(c_id[k])

			    # convert voxel-# to pixel scale
			    scaled_x = int(x/self.image_scaling) + self.img_min[0]
			    scaled_z = int(z/self.image_scaling) + self.img_min[2]
			    scaled_y = int(y/self.image_scaling) + self.img_min[1]

			    pos = (scaled_x, scaled_y, scaled_z)

			    # get new-position
			    #pos = self.get_position(pos)
			    #print pos, '-->', self.get_position(pos)

			    #if (scaled_x < self.img_width  and scaled_x > 0 and
			    if (scaled_z < self.img_width  and scaled_z > 0 and
			    	scaled_y < self.img_height and scaled_y > 0 ) :

				p0.append(pos)

		    #####
		    focal_x = self.settings.camera_focal_point[0]
		    ix = (self.img_max[0] - self.img_min[0])*focal_x

		    for iz in range(self.img_min[2], self.img_max[2]) :
		    #for ix in range(self.img_min[0], self.img_max[0]) :
		    	for iy in range(self.img_min[1], self.img_max[1]) :

			    pixel = (ix, iz, iy)
			    RGB = self._get_Signal(p0, pixel)

			    tirfm_image.putpixel((iz, self.img_height-1-iy), RGB)


		    tirfm_image.save(image_file_name)



        def make_movie(self) :
            """
            Make a movie by FFmpeg
            Requirement : Install FFmpeg (http://ffmpeg.org/)
            """
            input_image_filename = os.path.join(self.settings.movie_image_file_dir, self.settings.movie_image_file_name_format)
    
            # Set FFMPEG command
            ffmpeg_command  = 'ffmpeg -sameq -r "%s"' % (str(int(self.settings.camera_fps)))
            ffmpeg_command += ' -y -i "%s/%s" ' % (self.settings.movie_image_file_dir, self.settings.movie_image_file_name_format)
            ffmpeg_command += self.settings.movie_filename
            #ffmpeg_command += ('" -vcodec rawvideo -pix_fmt yuyv422 ' + self._movie_filename)
    
            os.system(ffmpeg_command)



	def output_movie(self, num_div=1):
	    """
	    Output movie
	    """
	    self.output_frames(num_div=num_div)
	    self.make_movie()



	def get_plots(self, plot_filename=None, wave_length=600) :

    	    from matplotlib.backends.backend_pdf import PdfPages
    
    	    pp = PdfPages(plot_filename)
    
    	    #####
            fig_spec_wl = pylab.figure()
    
    	    # fluorophore excitation and emission
    	    pylab.fill(self.settings.fluorophore_wavelength, self.settings.fluoex_eff, color='lightblue', label='Fluorophore Ex')
    	    pylab.fill(self.settings.fluorophore_wavelength, self.settings.fluoem_eff, color='pink', label='Fluorophore Em')
    
    	    # excitation filter
    	    if self.settings.excitation_switch == True :
    
    	        pylab.plot(self.settings.fluorophore_wavelength, self.settings.excitation_eff, color='blue', label='Excitation Filter', linewidth=2)
    
    	    # dichroic mirror
    	    if self.settings.dichroic_switch == True :
    
    	        pylab.plot(self.settings.fluorophore_wavelength, self.settings.dichroic_eff, color='green', label='Dichroic Mirror', linewidth=2)
    
            # emission filter
    	    if self.settings.emission_switch == True :
    
    	        pylab.plot(self.settings.fluorophore_wavelength, self.settings.emission_eff, color='red', label='Emission Filter', linewidth=2)
    
    	    pylab.axis([300, 800, 0, 100])
            pylab.xlabel('Wave Length [nm]')
            pylab.ylabel('Transmission Efficiency [%]')
            pp.savefig(fig_spec_wl)
    
    	    ######
    	    fig_spec_pos = pylab.figure()
    	    pylab.plot(self.settings.fluorophore_radial, self.settings.fluorophore_psf[0], color='pink', label='Radial PSF', linewidth=2)
    	    pylab.xlabel('Radial Position [nm]')
    	    pylab.ylabel('Photon [#]')
    	    pylab.title('Fluorophore : ' + self.settings.fluorophore_type)
    	    pp.savefig(fig_spec_pos)

            ######
            fig_spec_cont = pylab.figure()

            spec_scale = numpy.linspace(0, self.settings.fluorophore_psf[0][0], 101, endpoint=True)
            X, Y = numpy.meshgrid(self.settings.fluorophore_radial, self.settings.fluorophore_depth)

            pylab.contour(X, Y,  self.settings.fluorophore_psf, spec_scale, linewidth=0.1, color='k')
            pylab.contourf(X, Y, self.settings.fluorophore_psf, spec_scale, cmap=pylab.cm.jet)

            pylab.axis([0, 600, 0, 1000])
            pylab.xlabel('Radial [nm]')
            pylab.ylabel('Depth [nm]')
            pylab.title('Photon [#]')
            pp.savefig(fig_spec_cont)
    
    	    ######
    	    fig_qeff = pylab.figure()
	    pylab.plot(self.settings.fluorophore_wavelength, self.settings.camera_red,   color='red',   label='QE (Red)',   linewidth=2)
	    pylab.plot(self.settings.fluorophore_wavelength, self.settings.camera_green, color='green', label='QE (Green)', linewidth=2)
	    pylab.plot(self.settings.fluorophore_wavelength, self.settings.camera_blue,  color='blue',  label='QE (Blue)',  linewidth=2)
            pylab.axis([400, 1000, 0, 1.10])
            pylab.xlabel('Wave Length [nm]')
            pylab.ylabel('Quantum Efficiency')
            pylab.title('Camera : ' + self.settings.camera_type)
    	    pp.savefig(fig_qeff)

#            ###### SNR
#            self.settings.set_SNR()
#    
#    	    ###### SNR (wave_length)
#    	    index = int(wave_length) - self.settings.fluorophore_wavelength[0]
#    
#            fig_snr_wl= pylab.figure()
#            pylab.loglog(self.settings.photon_number, self.settings.absolute_snr[index])
#            pylab.plot(self.settings.photon_number, self.settings.ideal_snr, color='purple', label='Perfect', linewidth=2)
#            pylab.plot(self.settings.photon_number, self.settings.absolute_snr[index], color='red', label='SNR', linewidth=2)
#    
#            pylab.axis([self.settings.photon_number[0], self.settings.photon_number[-1], 0.01, self.settings.ideal_snr[-1]])
#            pylab.xlabel('Input Signal [photons/pixel]')
#            pylab.ylabel('SNR')
#            pylab.title('%d nm' % (int(wave_length)))
#            pp.savefig(fig_snr_wl)
#    
#            ###### Relative SNR (wave_length)
#            fig_rsnr_wl= pylab.figure()
#            pylab.semilogx(self.settings.photon_number)
#            pylab.plot(self.settings.photon_number, self.settings.ideal_relsnr, color='purple', label='Perfect', linewidth=2)
#            pylab.plot(self.settings.photon_number, self.settings.relative_snr[index], color='red', label='SNR', linewidth=2)
#    
#            pylab.axis([1, self.settings.photon_number[-1], 0, 1.10])
#            pylab.xlabel('Input Signal [photons/pixel]')
#            pylab.ylabel('Relative SNR')
#            pylab.title('%d nm' % (wave_length))
#            pp.savefig(fig_rsnr_wl)
#    
#            ###### SNR (Contour)
#            fig_snr = pylab.figure()
#            pylab.semilogx(self.settings.photon_number)
#    
#            snr_scale = numpy.linspace(0, self.settings.ideal_snr[-1], 21, endpoint=True)
#    	    X, Y = numpy.meshgrid(self.settings.photon_number, self.settings.fluorophore_wavelength)
#    
#            pylab.contour(X, Y, self.settings.absolute_snr,  snr_scale, linewidth=0.1, color='k')
#            pylab.contourf(X, Y, self.settings.absolute_snr, snr_scale, cmap=pylab.cm.jet)
#    
#    	    pylab.colorbar(ticks=snr_scale)
#            pylab.axis([self.settings.photon_number[0], self.settings.photon_number[-1], \
#    		self.settings.fluorophore_wavelength[100], self.settings.fluorophore_wavelength[800-self.settings.fluorophore_wavelength[0]]])
#            pylab.xlabel('Input Signal [photons/pixel]')
#            pylab.ylabel('Wave length [nm]')
#            pylab.title('SNR')
#            pp.savefig(fig_snr)
#    
#            ###### Relative SNR (Contour)
#            fig_rsnr = pylab.figure()
#            pylab.semilogx(self.settings.photon_number)
#    
#            rsnr_scale = numpy.linspace(0, 1, 21, endpoint=True)
#            X, Y = numpy.meshgrid(self.settings.photon_number, self.settings.fluorophore_wavelength)
#            pylab.contour(X, Y, self.settings.relative_snr,  rsnr_scale, linewidth=0.1, color='k')
#            pylab.contourf(X, Y, self.settings.relative_snr, rsnr_scale, cmap=pylab.cm.jet)
#    
#    	    pylab.colorbar(ticks=rsnr_scale)
#            pylab.axis([self.settings.photon_number[0], self.settings.photon_number[-1], \
#                    self.settings.fluorophore_wavelength[100], self.settings.fluorophore_wavelength[800-self.settings.fluorophore_wavelength[0]]])
#            pylab.xlabel('Input Signal [photons/pixel]')
#            pylab.ylabel('Wave length [nm]')
#            pylab.title('Relative SNR')
#            pp.savefig(fig_rsnr)
    
    	    pp.close()
    
    
###################################################################
#
#	Rotational Matrix
#
###################################################################
def trans_mat(v):
    """
    create rotation matrix for transform arbitrary vector
    to z-axis.
    """
    # rot_y on x-z plane
    ty = numpy.sign(v[0])*numpy.arccos( v[2]/numpy.sqrt(v[0]**2+v[2]**2) )
    ry = rot_y(ty)
    vym = numpy.dot(ry,v)
    vy = vym.A[0]
    
    # rot_x on y-z plane
    tx = -numpy.sign(vy[1])*numpy.arccos( vy[2]/numpy.sqrt(vy[1]**2+vy[2]**2) )
    rx = rot_x(tx)
    
    return numpy.dot(rx,ry)


def rot_x(t):
    """
    rotation matrix on x-axis
    """
    rx = numpy.matrix([
    [1.0, 0.0, 0.0],
    [0.0, numpy.cos(t), numpy.sin(t)],
    [0.0,-numpy.sin(t), numpy.cos(t)]
    ])

    return rx


def rot_y(t):
    """
    rotation matrix on y-axis
    """
    ry = numpy.matrix([
    [ numpy.cos(t), 0.0,-numpy.sin(t)],
    [ 0.0, 1.0, 0.0],
    [ numpy.sin(t), 0.0, numpy.cos(t)]
    ])

    return ry


def rot_z(t):
    """
    rotation matrix on z-axis
    """
    rz = numpy.matrix([
    [ numpy.cos(t), numpy.sin(t), 0.0],
    [-numpy.sin(t), numpy.cos(t), 0.0],
    [ 0.0, 0.0, 1.0]
    ])

    return rz


