"""

epifm_handler.py

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

from scipy.special import j0
from scipy.misc    import toimage

from matplotlib.backends.backend_pdf import PdfPages

IMAGE_SIZE_LIMIT=3000


class VisualizerError(Exception):

    "Exception class for visualizer"

    def __init__(self, info):
        self.__info = info

    def __repr__(self):
        return self.__info

    def __str__(self):
        return self.__info



class EPIFMConfigs() :

    '''
    EPIFM Visualization setting class

	Light Source
	Beam Expander
	Fluorophore
	Mirror
	Dichroic Mirror
	Excitation Filter
	Emission Filter
	Tube Lenz
	Detector
    '''

    def __init__(self, user_configs_dict = None):

        # default setting
        configs_dict = parameter_configs.__dict__.copy()

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


    def _set_data(self, key, val) :
    
        if val != None:
            setattr(self, key, val)



    def set_Particles(self, p_list = None) :

	print '--- Particles List :'

        self._set_data('particles_list', p_list)

	print p_list




    def set_LightSource(self,  source_type = None,
				wave_mode = None,
				M2_factor = None,
				wave_length = None,
                                power = None,
                                radius = None ) :

        self._set_data('source_switch', True)
        self._set_data('source_type', source_type)
        self._set_data('source_wavemode', wave_mode)
        self._set_data('source_M2factor', M2_factor)
        self._set_data('source_wavelength', wave_length)
        self._set_data('source_power', power)
        self._set_data('source_radius', radius)

        print '--- Light Source :', self.source_type
        print '\tWave Mode : ', self.source_wavemode
        print '\tM2 factor = ', self.source_M2factor
        print '\tWave Length = ', self.source_wavelength, 'nm'
        print '\tTotal Power = ', self.source_power, 'W (=joule/sec)'
        print '\t1/e2 Radius = ', self.source_radius, 'm'




    def set_BeamExpander(self, expander_type = None,
                                pinhole_radius = None,
                                focal_length1 = None,
                                focal_length2 = None) :

        self._set_data('expander_type', expander_type)
        self._set_data('expander_pinhole_radius', pinhole_radius)
        self._set_data('expander_focal_length1', focal_length1)
        self._set_data('expander_focal_length2', focal_length2)

        print '--- Beam Expander :', self.expander_type
        print '\tPinhole Radius = ', self.expander_pinhole_radius, 'm'
        print '\tFocal Length 1 = ', self.expander_focal_length1, 'm'
        print '\tFocal Length 2 = ', self.expander_focal_length2, 'm'




    def set_ExcitationFilter(self, excitation = None) :

        print '--- Excitation Filter :'

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

	    index = self.psf_wavelength - self.wave_length[0]
	    self.fluoex_eff[index] = 100
	    self.fluoem_eff[index] = 100

            print '\tWave Length   = ', self.psf_wavelength, 'nm'
            print '\tLateral Width = ', self.psf_width[0], 'nm'
            print '\tAxial Width = ', self.psf_width[1], 'nm'
	    print '\tLateral Cutoff = ', self.psf_cutoff[0], 'nm'
            print '\tAxial Cutoff = ', self.psf_cutoff[1], 'nm'



	elif (fluorophore_type == 'Point-like') :

            print '--- Fluorophore : Point Spreading Function [%s]' % (fluorophore_type)

            self._set_data('fluorophore_type', fluorophore_type)
            self._set_data('psf_wavelength', wave_length)
            self._set_data('psf_file_name_format', file_name_format)

            index = self.psf_wavelength - self.wave_length[0]
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
	    self.fluoex_eff = self.set_efficiency(fluorophore_excitation)
	    self.fluoem_eff = self.set_efficiency(fluorophore_emission)

	    index_ex = self.fluoex_eff.index(max(self.fluoex_eff))
	    index_em = self.fluoem_eff.index(max(self.fluoem_eff))

	    #### for temporary
            self._set_data('fluorophore_type', fluorophore_type)
	    self._set_data('psf_wavelength', self.wave_length[index_em])
            self._set_data('psf_file_name_format', file_name_format)

            self.fluoex_eff[index_ex] = 100
            self.fluoem_eff[index_em] = 100

            print '\tExcitation : Wave Length   = ', self.wave_length[index_ex], 'nm'
            print '\tEmission   : Wave Length   = ', self.psf_wavelength, 'nm'


	# Normalization
	norm = sum(self.fluoex_eff)
	self.fluoex_norm = numpy.array(self.fluoex_eff)/norm

	norm = sum(self.fluoem_eff)
	self.fluoem_norm = numpy.array(self.fluoem_eff)/norm



    def set_Objective(self, 
			NA = None,
			Ng = None,
			Nm = None,
			focal_length = None,
			efficiency = None,
			thickness  = None
			) :

	print '--- Objective :'

        self._set_data('objective_switch', True)
        self._set_data('objective_NA', NA)
        self._set_data('objective_Ng', Ng)
        self._set_data('objective_Nm', Nm)
        self._set_data('objective_focal_length', focal_length)
        self._set_data('objective_efficiency', efficiency)
        self._set_data('objective_glassthickness', thickness)
	self._set_angle()

	print '\tNA = ', self.objective_NA
	print '\tN(glass)  = ', self.objective_Ng
	print '\tN(medium) = ', self.objective_Nm
	print '\tFocal Length  = ', self.objective_focal_length, 'm'
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



    def set_ScanLens(self, focal_length = None) :

        print '--- Scan Lens :'

        self._set_data('scanlens_switch', True)
        self._set_data('scanlens_focal_length', focal_length)

        print '\tFocal Length = ', self.scanlens_focal_length, 'm'



    def set_TubeLens(self, focal_length = None) :

        print '--- Tube Lens :'
	
        self._set_data('tubelens_switch', True)
        self._set_data('tubelens_focal_length', focal_length)

        print '\tFocal Length = ', self.tubelens_focal_length, 'm'




    def set_Detector(self, detector = None,
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

        print '--- Detector :'

        
        try:
	    filename = './catalog/detector/' + detector + '.csv'

            csvfile = open(filename)
            lines = csvfile.readlines()

            header = lines[0:12]
            data   = lines[13:]

            detector_header = []
            detector_QEdata = []

            for i in range(len(header)) :
                dummy  = header[i].split('\r\n')
                a_data = dummy[0].split(',')
                detector_header.append(a_data)
                #print '\t', a_data

	    image_size	 = (int(detector_header[5][1]), int(detector_header[5][1]))
	    pixel_length = float(detector_header[6][1])
	    sat_charge   = float(detector_header[7][1])
	    readout	 = float(detector_header[8][1])
	    dark_current = float(detector_header[9][1])
	    excess	 = float(detector_header[10][1])
	    emgain	 = float(detector_header[11][1])

            for i in range(len(data)) :
                dummy0 = data[i].split('\r\n')
                a_data = dummy0[0].split(',')
                detector_QEdata.append(a_data)

	    self.detector_qeff  = self.set_efficiency(detector_QEdata)
	    #self.detector_blue  = self.set_efficiency(detector_QEdata, 1)
	    #self.detector_green = self.set_efficiency(detector_QEdata, 2)
	    #self.detector_red   = self.set_efficiency(detector_QEdata, 3)



        except Exception:

	    print 'Error : ', filename, ' file is NOT found'
	    print '\tUse Perfect detector (Noise free)\n'

	    detector = 'Perfect'
            #exit()

        self._set_data('detector_switch', True)
        self._set_data('detector_type', detector)
        self._set_data('detector_image_size', image_size)
        self._set_data('detector_pixel_length', pixel_length)
        self._set_data('detector_focal_point', focal_point)
        self._set_data('detector_base_position', base_position)
        self._set_data('detector_zoom', zoom)
        self._set_data('detector_start_time', start_time)
        self._set_data('detector_end_time', end_time)
        self._set_data('detector_fps', fps)
        self._set_data('detector_exposure_time', exposure_time)
        self._set_data('detector_sat_charge', sat_charge)
        self._set_data('detector_ADC_bit', ADC_bit)
        self._set_data('detector_ADC_const', ADC_const)
        self._set_data('detector_ADC_offset', ADC_offset)
        self._set_data('detector_dark_current', dark_current)
        self._set_data('detector_readout', readout)
        self._set_data('detector_excess', excess)
        self._set_data('detector_emgain', emgain)

	print '\tDetector Type : ', self.detector_type
        print '\tImage Size  = ', self.detector_image_size[0], 'x', self.detector_image_size[1]
        print '\tPixel Size  = ', self.detector_pixel_length, 'm/pixel'
        print '\tFocal Point = ', self.detector_focal_point
        print '\tPosition    = ', self.detector_base_position
        print '\tZoom	     = ', self.detector_zoom
        print '\tStart Time  = ', self.detector_start_time, 'sec'
        print '\tEnd   Time  = ', self.detector_end_time, 'sec'
        print '\tFrame Rate  = ', self.detector_fps, 'frames/sec'
        print '\tExposure Time = ', self.detector_exposure_time, 'sec'
        print '\tSat. Charge   = ', self.detector_sat_charge, 'electron'
        print '\tA/D Converter = ', self.detector_ADC_bit, 'bit'
        print '\tADC Const  = ', self.detector_ADC_const, 'electron/ADC count'
        print '\tADC Offset = ', self.detector_ADC_offset, 'ADC count'
        print '\tDark Current = ', self.detector_dark_current
        print '\tReadout = ', self.detector_readout
        print '\tExcess	 = ', self.detector_excess
        print '\tEM gain = ', self.detector_emgain



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

        N = len(self.wave_length)
        #efficiency = numpy.array([0.0 for i in range(N)])
        efficiency = [0.0 for i in range(N)]

        for i in range(N) :
            wl = self.wave_length[i]

            for j in range(len(array)) :

                length = float(array[j][0])
                eff = float(array[j][index])

                if (length/wl == 1) :
                    efficiency[i] = eff


        return efficiency




    def set_Optical_path(self) :

	# (1) Illumination path : Light source --> Sample
	self.set_Illumination_path()

	# (2) Scattering matrix : photon-molecule interaction
	self.set_SMatrix()

	# (3) Detection path : Sample --> Detector
	self.set_Detection_path()



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
	penetration_depth = float('inf')

        psf_r  = numpy.array(map(lambda x : 1.00, r))
        psf_z  = numpy.exp(-z/penetration_depth)
        psf_pd = numpy.array(map(lambda x : psf_r*x, psf_z))

        # Illumination PSF : Total
        self.source_psf = F0*bsf*psf_pd



    def set_SMatrix(self) :

        #################################################################
        #
        # Note : Scattering Matrix of photon-matter interactions
        #
        #       Using random number generator to simulate
        #       the energy transition of molecular vibrational states
        #
        #################################################################

        r = self.radial
        z = self.depth

        wave_length = self.wave_length

	M = 1.00e-6

	self.scattering_matrix = M

        #################################################################


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

        view = self.detector_pixel_length/(2.0*voxel_radius)
        zoom = self.detector_zoom

        self.image_scaling = view/(Mag*zoom)

        # Detector PSF
        self.set_PSF_detector(Mag)





    def set_PSF_detector(self, Mag) :

        r = self.radial
        z = self.depth

        wave_length = self.psf_wavelength

        # Fluorophores Emission Intensity (wave_length)
        I = self.fluoem_norm

        # Photon Transmission Efficiency
        if (self.dichroic_switch == True) :
            I = I*0.01*self.dichroic_eff

        if (self.emission_switch == True) :
            I = I*0.01*self.emission_eff

        # Detector : Quantum Efficiency
	# Convertion of photons to photoelectrons
        I = I*self.detector_qeff

	# For normalization
	norm = map(lambda x : True if x > 1e-4 else False, I)

        # PSF : Focal Depth
        n1 = self.objective_Ng
        n2 = self.objective_Nm

        depth_1st = n1/(2.0*n2**2)*wave_length
        depth_2nd = n1/(7*Mag*n2)

        focal_depth = depth_1st + depth_2nd

	psf_r  = numpy.array(map(lambda x : 1.00, r))
        psf_z  = numpy.exp(-0.5*(z/focal_depth)**2)

	psf_fd = numpy.array(map(lambda x : psf_r*x, psf_z))


	# PSF : Fluorophore
	psf_fl = None

	if (self.fluorophore_type == 'Gaussian') :

            I0 = 1.0
            Ir = sum(map(lambda x : x*numpy.exp(-0.5*(r/self.psf_width[0])**2), norm))
            Iz = sum(map(lambda x : x*numpy.exp(-0.5*(z/self.psf_width[1])**2), norm))

	    psf_fl = numpy.array(map(lambda x : I0*Ir*x, Iz))


        elif (self.fluorophore_type == 'Point-like') :

            I0 = 1.0
            Ir = sum(map(lambda x : x*numpy.array(map(lambda x : 1.00 if x == 0 else 0.00, r)), norm))
            Iz = sum(map(lambda x : x*numpy.array(map(lambda x : 1.00 if x == 0 else 0.00, z)), norm))

            psf_fl = numpy.array(map(lambda x : I0*Ir*x, Iz))


	else :

	    # make the norm and wave_length array shorter
#	    re_wavelength = []
#
#	    for i in range(len(norm)) :
#
#		if norm[i] is True :
#		    re_wavelength.append(wave_length[i])
#
#	    psf = self.get_PSF_fluorophore(r, z, numpy.array(re_wavelength))
	    psf_fl = self.get_PSF_fluorophore(r, z, wave_length)

	self.fluorophore_psf = psf_fl*psf_fd



    def get_PSF_fluorophore(self, r, z, wave_length) :

	# set Numerical Appature
	NA = self.objective_NA
	N = 80

	# set alpha and gamma consts
	k = 2.0*numpy.pi/wave_length
	alpha = k*NA
	gamma = k*(NA/2)**2

	# set rho parameters
        drho = 1.0/N
        rho = numpy.array([(i+1)*drho for i in range(N)])

        J0 = numpy.array(map(lambda x : j0(x*alpha*rho), r))
        Y  = numpy.array(map(lambda x : 2*numpy.exp(-2*1.j*x*gamma*rho**2)*rho*drho, z))
#	J0 = numpy.array(map(lambda y : map(lambda x : j0(x*y*rho), r), alpha))
#	Y  = numpy.array(map(lambda y : map(lambda x : 2*numpy.exp(-2*1.j*x*y*rho**2)*rho*drho, z), gamma))

        I  = numpy.array(map(lambda x : x*J0, Y))
        I_sum = I.sum(axis=2)

        psf = numpy.array(map(lambda x : abs(x)**2, I_sum))

#	for i in range(len(wave_length)) :
#
#	    I  = numpy.array(map(lambda x : x*J0[i], Y[i]))
#	    I_sum = I.sum(axis=2)
#	    I_abs = map(lambda x : abs(x)**2, I_sum)
#
#	    if (i > 0) : psf += I_abs
#	    else : psf = I_abs

	return psf



#    def set_PSF_detector(self, scaling_factor) :
#
#        # Detector : Quantum Efficiency
#        self.fluorophore_rgb[:,2] = map(lambda x : sum(x), map(lambda x : x*self.detector_blue,  N_ph))
#        self.fluorophore_rgb[:,1] = map(lambda x : sum(x), map(lambda x : x*self.detector_green, N_ph))
#        self.fluorophore_rgb[:,0] = map(lambda x : sum(x), map(lambda x : x*self.detector_red,   N_ph))
#
#        # count photoelectrons
#        N_b = sum(map(lambda x : x*self.detector_blue,  N_ph))
#        N_g = sum(map(lambda x : x*self.detector_green, N_ph))
#        N_r = sum(map(lambda x : x*self.detector_red,   N_ph))
#
#        voxel_radius = self.spatiocyte_VoxelRadius
#        pixel_length = int(2.0*voxel_radius*scaling_factor/1e-9)
#
#	Nr = len(self.radial)/pixel_length
#	#Nz = len(self.depth)/pixel_length
#
#	psf_in_pixel = numpy.array([[0.00 for i in range(Nr)] for j in range(len(self.depth))])
#
#	for i in range(len(self.depth)) :
#
#	    for j in range(pixel_length) :
#
#		array_r = self.fluorophore_psf[i,:Nr*pixel_length].reshape((Nr, pixel_length))
#		psf_in_pixel[i] = map(lambda x : sum(x), array_r)
#
#	psf_in_pixel = psf_in_pixel/psf_in_pixel[0][0]



    def get_Noise(self, signal) :

	Nr = self.detector_readout
	DC = self.detector_dark_current
	Fn = numpy.sqrt(self.detector_excess)
	M  = self.detector_emgain

	# Noise calculation defined by HAMAMATSU Photonics
	sigma2 = Fn**2*(signal + DC*self.detector_exposure_time)+ (Nr/M)**2
	noise  = numpy.sqrt(sigma2)

	return noise



    def set_SNR(self) :

	for j in range(len(self.wave_length)) :

            #QEff = self.detector_green[j]
            QEff = self.detector_qeff[j]

	    # get SNR and relative SNR
	    self.absolute_snr[j] = (QEff*self.photon_number)/map(lambda x : self.get_Noise(QEff*x), self.photon_number)
	    self.relative_snr[j] = map(lambda x, y : x/y, self.absolute_snr[j], self.ideal_snr)



class EPIFMVisualizer() :

	'''
	EPIFM Visualization class of e-cell simulator
	'''

	def __init__(self, configs=EPIFMConfigs()) :

		assert isinstance(configs, EPIFMConfigs)
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
                set Image Size and Boundary
                """
                self.img_width  = int(self.configs.detector_image_size[0])
                self.img_height = int(self.configs.detector_image_size[1])

                if self.img_width > IMAGE_SIZE_LIMIT or self.img_height > IMAGE_SIZE_LIMIT :
                        raise VisualizerErrror('Image size is bigger than the limit size')

                """
                set Optical path from light source to detector
                """
		self.configs.set_Optical_path()



        def __del__(self):

            if self.configs.movie_cleanup_image_file_dir :

                for parent_dir, dirs, files in os.walk(self.configs.movie_image_file_dir, False) :
                    for file in files :
                        os.remove(os.path.join(parent_dir, file))

                    os.rmdir(parent_dir)



	def get_coordinate(self, aCoord) :

	        """
		get (column, layer, row) coordinate
	        """
		start_coord = self.configs.spatiocyte_theStartCoord
		row_size    = self.configs.spatiocyte_theRowSize
		layer_size  = self.configs.spatiocyte_theLayerSize
		col_size    = self.configs.spatiocyte_theColSize

                aGlobalCol   =  (aCoord-start_coord)/(row_size*layer_size)
                aGlobalLayer = ((aCoord-start_coord)%(row_size*layer_size))/row_size
                aGlobalRow   = ((aCoord-start_coord)%(row_size*layer_size))%row_size

                """
		get (x, y, z) coordinate
                """
		norm_voxel_radius = self.configs.spatiocyte_theNormalizedVoxelRadius

                theHCPk = norm_voxel_radius/math.sqrt(3.0)
                theHCPh = norm_voxel_radius*math.sqrt(8.0/3.0)
                theHCPl = norm_voxel_radius*math.sqrt(3.0)

	        point_y = (aGlobalCol%2)*theHCPk + theHCPl*aGlobalLayer
	        point_z = aGlobalRow*2*norm_voxel_radius + ((aGlobalLayer+aGlobalCol)%2)*norm_voxel_radius
	        point_x = aGlobalCol*theHCPh

		return point_x, point_y, point_z



	def get_position(self, pos) :

                # normal vector
                norm  = map(operator.sub, self.configs.detector_base_position, self.configs.detector_focal_point)
                norm_len2 = norm[0]**2 + norm[1]**2 + norm[2]**2
                norm_len = math.sqrt(norm_len2)

                # detector's focal position
                focal = numpy.array(self.configs.detector_focal_point)

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




        def get_signal(self, p_i, p_b, p_0) :

		# set focal point
                x_0, y_0, z_0 = p_0

                # set source center
                x_b, y_b, z_b = p_b

		# set particle position
                x_i, y_i, z_i = p_i

		#####
		r = self.configs.radial
		d = self.configs.depth

                z = numpy.linspace(-r[-1], +r[-1], 2*len(r)-1)
                y = numpy.linspace(-r[-1], +r[-1], 2*len(r)-1)

                Z, Y = numpy.meshgrid(z, y)
		R = numpy.sqrt(Z**2 + Y**2)


		# get PSF (Light Source)
		dd = abs(x_b - x_0)

                if (dd < len(d)) :
		    source_depth = dd
                else :
		    source_depth = d[-1]

		# beam lateral position
		rr = numpy.sqrt((y_i-y_0)**2 + (z_i-z_0)**2)

		if (rr < len(r)) :
		    source_radius = rr
		else :
		    source_radius = r[-1]

		source_psf = self.configs.source_psf[int(source_depth)][int(source_radius)]


		# get scattering matrix (y-z)
		SMatrix = self.configs.scattering_matrix


                # get PSF (Fluorophore)
		if (abs(x_i - x_0) >= len(d)) :
		    fluo_depth = d[-1]
		else :
		    fluo_depth = abs(x_i - x_0)

                fluo_psf = copy.copy(R)
                fluo_array = self.configs.fluorophore_psf[int(fluo_depth)]

		for i in range (len(R)) :
		    for j in range(len(R[i])) :

			rr = R[i][j]

                        if (rr < len(fluo_array)) :
			    fluo_psf[i][j] = fluo_array[int(rr)]
                        else :
			    fluo_psf[i][j] = 0


                # signal conversion : PSF out = PSF(Fluorophore) * SMatrix * PSF(Beam)
		signal = fluo_psf * SMatrix * source_psf


                return signal



	def overwrite(self, plane, signal, p_i) :

                # particle position
                x_i, y_i, z_i = p_i

		# z-axis
		Nz_plane  = len(plane)
		Nz_signal = len(signal)
		Nr = len(self.configs.radial)

                z_to   = z_i + Nr
                z_from = z_i - Nr

                if (z_to > Nz_plane) :

                    dz_to = z_to - Nz_plane

                    z0_to = int(Nz_plane)
                    zi_to = int(Nz_signal - dz_to)

                else :

                    dz_to = Nz_plane - (z_i + Nr)

                    z0_to = int(Nz_plane - dz_to)
                    zi_to = int(Nz_signal)

                if (z_from < 0) :

                    dz_from = abs(z_from)

                    z0_from = 0
                    zi_from = int(dz_from)

                else :

                    dz_from = z_from

                    z0_from = int(dz_from)
                    zi_from = 0

                ddz = (z0_to - z0_from) - (zi_to - zi_from)

                if (ddz > 0) : z0_to = z0_to - ddz
                if (ddz < 0) : zi_to = zi_to - ddz

                # y-axis
                Ny_plane  = plane.size/Nz_plane
                Ny_signal = signal.size/Nz_signal

		y_to   = y_i + Nr
		y_from = y_i - Nr

                if (y_to > Ny_plane) :

                    dy_to = y_to - Ny_plane

                    y0_to = int(Ny_plane)
                    yi_to = int(Ny_signal - dy_to)

                else :

                    dy_to = Ny_plane - (y_i + Nr)

                    y0_to = int(Ny_plane - dy_to)
                    yi_to = int(Ny_signal)

                if (y_from < 0) :

                    dy_from = abs(y_from)

                    y0_from = 0
                    yi_from = int(dy_from)

                else :

                    dy_from = y_from

                    y0_from = int(dy_from)
                    yi_from = 0

		ddy = (y0_to - y0_from) - (yi_to - yi_from)

		if (ddy > 0) : y0_to = y0_to - ddy
		if (ddy < 0) : yi_to = yi_to - ddy

		# add to plane
                plane[z0_from:z0_to, y0_from:y0_to] += signal[zi_from:zi_to, yi_from:yi_to]


		return plane


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


                    #fig = pylab.figure()
                    #spec_scale = numpy.linspace(numpy.amin(cell), numpy.amax(cell), 200, endpoint=True)
                    #pylab.contour(Z, Y,  cell, spec_scale, linewidth=0.1, color='k')
                    #pylab.contourf(Z, Y, cell, cmap=pylab.cm.jet)
                    #pylab.show()
                    #exit()

                    self.detector_output(cell)

		    time  += frame_interval
		    count += 1



	def detector_output(self, cell) :

		# Detector Output
                voxel_size = (2.0*self.configs.spatiocyte_VoxelRadius)/1e-9

		Nw_pixel = self.img_width
		Nh_pixel = self.img_height

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

		for i in range(Nh_cell/Np) :
		    for j in range(Nw_cell/Np) :

			signal = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])
			cell_pixel[i][j] = self.A2D_converter(signal)


                # image in pixel-scale
		camera_pixel = numpy.zeros(shape=(Nw_pixel, Nh_pixel))


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

                # save image to file
                max_bit = 2**self.configs.detector_ADC_bit

                #toimage(camera_pixel, cmin=0, cmax=max_bit).save(self.image_file_name)
                toimage(camera_pixel, cmin=0, cmax=numpy.amax(camera_pixel)).save(self.image_file_name)

		#z = numpy.linspace(0, Nw_pixel-1, Nw_pixel)
		#y = numpy.linspace(0, Nh_pixel-1, Nh_pixel)
		#Z, Y = numpy.meshgrid(z, y)

		#fig = pylab.figure()
		#spec_scale = numpy.linspace(numpy.amin(camera_pixel), numpy.amax(camera_pixel), 200, endpoint=True)
		#pylab.contour(Z, Y,  camera_pixel, spec_scale, linewidth=0.1, color='k')
		#pylab.contourf(Z, Y, camera_pixel, cmap=pylab.cm.jet)
		#pylab.show()
		#exit()



	def A2D_converter(self, signal) :

            # get noise
            noise = self.configs.get_Noise(signal)

            # convert photoelectron to ADC counts (Grayscale)
            k = self.configs.detector_ADC_const
            ADC0 = self.configs.detector_ADC_offset
	    Gain = self.configs.detector_emgain

            ADC  = signal
            #ADC  = Gain*numpy.random.poisson(signal, None)
            #ADC  = Gain*numpy.random.poisson((signal+noise)/k + ADC0, None)

	    return ADC



        def make_movie(self) :
            """
            Make a movie by FFmpeg
            Requirement : Install FFmpeg (http://ffmpeg.org/)
            """
            input_image_filename = os.path.join(self.configs.movie_image_file_dir, self.configs.movie_image_file_name_format)
    
            # Set FFMPEG command
            ffmpeg_command  = 'ffmpeg -sameq -r "%s"' % (str(int(self.configs.detector_fps)))
            ffmpeg_command += ' -y -i "%s/%s" ' % (self.configs.movie_image_file_dir, self.configs.movie_image_file_name_format)
            ffmpeg_command += self.configs.movie_filename
            #ffmpeg_command += ('" -vcodec rawvideo -pix_fmt yuyv422 ' + self._movie_filename)
    
            os.system(ffmpeg_command)



	def output_movie(self, num_div=1):
	    """
	    Output movie
	    """
	    self.output_frames(num_div=num_div)
	    self.make_movie()



	def get_plots(self, plot_filename=None, wave_length=600) :

            pdf_page = PdfPages(plot_filename)

	    self.plot_LightSource(pdf_page)
	    self.plot_Fluorophore(pdf_page)
	    self.plot_Detector(pdf_page)
	    self.plot_SNR(pdf_page, wave_length)

	    pdf_page.close()



	def plot_LightSource(self, pdf_page) :

            #####
#            fig_spec_wl = pylab.figure()
#
#            # beam spectrum
#            pylab.fill(self.configs.wave_length, self.configs.fluoex_eff, color='lightblue', label='Incident Beam')
#
#            # excitation filter
#            if self.configs.excitation_switch == True :
#
#                pylab.plot(self.configs.wave_length, self.configs.excitation_eff, color='blue', label='Excitation Filter', linewidth=2)
#
#
#            pylab.axis([300, 800, 0, 100])
#            pylab.xlabel('Wave Length [nm]')
#            pylab.ylabel('Beam Intensity : %f [W/cm^2]' % (self.configs.source_intensity))
#
#            pdf_page.savefig(fig_spec_wl)

            ######
            fig_spec_pos = pylab.figure()

            pylab.plot(self.configs.radial, self.configs.source_psf[0], color='pink', label='Radial PSF', linewidth=2)
            pylab.xlabel('Radial Position [nm]')
            pylab.ylabel('Photon [#]')
            pylab.title(self.configs.source_type)

            pdf_page.savefig(fig_spec_pos)

            ######
            fig_spec_cont = pylab.figure()

            spec_scale = numpy.linspace(0, self.configs.source_psf[0][0], 101, endpoint=True)
            X, Y = numpy.meshgrid(self.configs.radial, self.configs.depth)

            pylab.contour(X, Y,  self.configs.source_psf, spec_scale, linewidth=0.1, color='k')
            pylab.contourf(X, Y, self.configs.source_psf, spec_scale, cmap=pylab.cm.jet)

            pylab.axis([0, self.configs.radial[-1], 0, self.configs.depth[-1]])
            pylab.xlabel('Radial [nm]')
            pylab.ylabel('Depth [nm]')
            pylab.title('Beam Photon [#]')

            pdf_page.savefig(fig_spec_cont)




	def plot_Fluorophore(self, pdf_page) :

    	    #####
            fig_spec_wl = pylab.figure()
    
    	    # fluorophore excitation and emission
    	    pylab.fill(self.configs.wave_length, self.configs.fluoex_eff, color='lightblue', label='Fluorophore Ex')
    	    pylab.fill(self.configs.wave_length, self.configs.fluoem_eff, color='pink', label='Fluorophore Em')
    
    	    # excitation filter
    	    if self.configs.excitation_switch == True :
    
    	        pylab.plot(self.configs.wave_length, self.configs.excitation_eff, color='blue', label='Excitation Filter', linewidth=2)
    
    	    # dichroic mirror
    	    if self.configs.dichroic_switch == True :
    
    	        pylab.plot(self.configs.wave_length, self.configs.dichroic_eff, color='green', label='Dichroic Mirror', linewidth=2)
    
            # emission filter
    	    if self.configs.emission_switch == True :
    
    	        pylab.plot(self.configs.wave_length, self.configs.emission_eff, color='red', label='Emission Filter', linewidth=2)
    
    	    pylab.axis([300, 800, 0, 100])
            pylab.xlabel('Wave Length [nm]')
            pylab.ylabel('Transmission Efficiency [%]')

            pdf_page.savefig(fig_spec_wl)
    
    	    ######
    	    fig_spec_pos = pylab.figure()

    	    pylab.plot(self.configs.radial, self.configs.fluorophore_psf[0], color='pink', label='Radial PSF', linewidth=2)
    	    pylab.xlabel('Radial Position [nm]')
    	    pylab.ylabel('Single Photon')
    	    pylab.title('Fluorophore : ' + self.configs.fluorophore_type)

    	    pdf_page.savefig(fig_spec_pos)

            ######
            fig_spec_cont = pylab.figure()

            spec_scale = numpy.linspace(0, self.configs.fluorophore_psf[0][0], 101, endpoint=True)
            X, Y = numpy.meshgrid(self.configs.radial, self.configs.depth)

            pylab.contour(X, Y,  self.configs.fluorophore_psf, spec_scale, linewidth=0.1, color='k')
            pylab.contourf(X, Y, self.configs.fluorophore_psf, spec_scale, cmap=pylab.cm.jet)

            pylab.axis([0, 600, 0, 1000])
            pylab.xlabel('Radial [nm]')
            pylab.ylabel('Depth [nm]')
            pylab.title('Photons')
            pdf_page.savefig(fig_spec_cont)




        def plot_Detector(self, pdf_page) :

    	    ######
    	    fig_qeff = pylab.figure()

	    #pylab.plot(self.configs.wave_length, self.configs.detector_red,   color='red',   label='QE (Red)',   linewidth=2)
	    #pylab.plot(self.configs.wave_length, self.configs.detector_green, color='green', label='QE (Green)', linewidth=2)
	    #pylab.plot(self.configs.wave_length, self.configs.detector_blue,  color='blue',  label='QE (Blue)',  linewidth=2)
	    pylab.plot(self.configs.wave_length, self.configs.detector_qeff,  color='black',  label='QE',  linewidth=2)
            pylab.axis([400, 1000, 0, 1.10])
            pylab.xlabel('Wave Length [nm]')
            pylab.ylabel('Quantum Efficiency')
	    pylab.title('Camera : ' + self.configs.detector_type)

    	    pdf_page.savefig(fig_qeff)




        def plot_SNR(self, pdf_page, wave_length=600) :

            ###### SNR
            self.configs.set_SNR()
    
    	    ###### SNR (wave_length)
    	    index = (numpy.abs(self.configs.wave_length - int(wave_length))).argmin()
    
            fig_snr_wl= pylab.figure()
            pylab.loglog(self.configs.photon_number, self.configs.absolute_snr[index])
            pylab.plot(self.configs.photon_number, self.configs.ideal_snr, color='purple', label='Perfect', linewidth=2)
            pylab.plot(self.configs.photon_number, self.configs.absolute_snr[index], color='red', label='SNR', linewidth=2)
    
            pylab.axis([self.configs.photon_number[0], self.configs.photon_number[-1], 0.01, self.configs.ideal_snr[-1]])
            pylab.xlabel('Input Signal [photons/pixel]')
            pylab.ylabel('SNR')
            pylab.title('%d nm' % (int(wave_length)))
            pdf_page.savefig(fig_snr_wl)
    
            ###### Relative SNR (wave_length)
            fig_rsnr_wl= pylab.figure()
            pylab.semilogx(self.configs.photon_number)
            pylab.plot(self.configs.photon_number, self.configs.ideal_relsnr, color='purple', label='Perfect', linewidth=2)
            pylab.plot(self.configs.photon_number, self.configs.relative_snr[index], color='red', label='SNR', linewidth=2)
    
            pylab.axis([1, self.configs.photon_number[-1], 0, 1.10])
            pylab.xlabel('Input Signal [photons/pixel]')
            pylab.ylabel('Relative SNR')
            pylab.title('%d nm' % (wave_length))
            pdf_page.savefig(fig_rsnr_wl)
    
            ###### SNR (Contour)
#            fig_snr = pylab.figure()
#            pylab.semilogx(self.configs.photon_number)
#    
#            snr_scale = numpy.linspace(0, self.configs.ideal_snr[-1], 21, endpoint=True)
#    	    X, Y = numpy.meshgrid(self.configs.photon_number, self.configs.wave_length)
#    
#            pylab.contour(X, Y, self.configs.absolute_snr,  snr_scale, linewidth=0.1, color='k')
#            pylab.contourf(X, Y, self.configs.absolute_snr, snr_scale, cmap=pylab.cm.jet)
#    
#    	    pylab.colorbar(ticks=snr_scale)
#            pylab.axis([self.configs.photon_number[0], self.configs.photon_number[-1], \
#    		self.configs.wave_length[0], self.configs.wave_length[-1]])
#            pylab.xlabel('Input Signal [photons/pixel]')
#            pylab.ylabel('Wave length [nm]')
#            pylab.title('SNR')
#            pdf_page.savefig(fig_snr)
#    
#            ###### Relative SNR (Contour)
#            fig_rsnr = pylab.figure()
#
#            pylab.semilogx(self.configs.photon_number)
#    
#            rsnr_scale = numpy.linspace(0, 1, 21, endpoint=True)
#            X, Y = numpy.meshgrid(self.configs.photon_number, self.configs.wave_length)
#            pylab.contour(X, Y, self.configs.relative_snr,  rsnr_scale, linewidth=0.1, color='k')
#            pylab.contourf(X, Y, self.configs.relative_snr, rsnr_scale, cmap=pylab.cm.jet)
#    
#    	    pylab.colorbar(ticks=rsnr_scale)
#            pylab.axis([self.configs.photon_number[0], self.configs.photon_number[-1], \
#                    self.configs.wave_length[0], self.configs.wave_length[-1]])
#            pylab.xlabel('Input Signal [photons/pixel]')
#            pylab.ylabel('Wave length [nm]')
#            pylab.title('Relative SNR')
#
#            pdf_page.savefig(fig_rsnr)


    
###################################################################
#
#	Rotational Matrix
#
###################################################################
def trans_mat(v):
    """
    create rotation matrix for transform arbitrary vector to z-axis.
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


