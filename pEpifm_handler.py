"""

epifm_handler.py

"""
import sys
import os
import copy
import tempfile
import math
import operator
import random
import h5py
import csv
import ctypes
import multiprocessing

#import pylab
import scipy
import numpy
import Image

import parameter_configs
import parameter_effects
from effects_handler import PhysicalEffects

from time import sleep

from scipy.special import j0
from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates
from scipy.misc    import toimage

#from matplotlib.backends.backend_pdf import PdfPages

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

	    index = (numpy.abs(self.wave_length - self.psf_wavelength)).argmin()

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
            self._set_data('psf_width', (10, 140))
            self._set_data('psf_file_name_format', file_name_format)

            index = (numpy.abs(self.wave_length - self.psf_wavelength)).argmin()

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



    def set_TubeLens1(self, focal_length = None) :

        print '--- Tube Lens :'
	
        self._set_data('tubelens_switch', True)
        self._set_data('tubelens_focal_length1', focal_length)

        print '\tFocal Length = ', self.tubelens_focal_length1, 'm'



    def set_TubeLens2(self, focal_length = None) :

        print '--- Tube Lens :'

        self._set_data('tubelens_switch', True)
        self._set_data('tubelens_focal_length2', focal_length)

        print '\tFocal Length = ', self.tubelens_focal_length2, 'm'




    def set_Detector(self, detector = None,
		   image_size = None,
                   pixel_length = None,
                   focal_point = None,
                   base_position = None,
                   zoom = None,
                   start_time = None,
                   end_time = None,
                   fps = None,
                   exposure_time = None,
		   ADC_satQ = None,
		   ADC_maxQ = None,
		   ADC_bit = None,
		   fpn_type = None,
		   fpn_offset = None,
		   fpn_noise  = None,
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

            header = lines[0:14]
            data   = lines[15:]

            detector_header = []
            detector_QEdata = []

            for i in range(len(header)) :
                dummy  = header[i].split('\r\n')
                a_data = dummy[0].split(',')
                detector_header.append(a_data)
                #print '\t', a_data

	    image_size	 = (int(detector_header[5][1]), int(detector_header[5][1]))
	    pixel_length = float(detector_header[6][1])
	    ADC_maxQ   = float(detector_header[7][1])
	    ADC_satQ   = float(detector_header[8][1])
	    ADC_bit 	 = float(detector_header[9][1])
	    fpn_offset	 = float(detector_header[10][1])
	    readout	 = float(detector_header[11][1])
	    dark_current = float(detector_header[12][1])
	    excess	 = float(detector_header[13][1])

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
        self._set_data('detector_ADC_maxQ', ADC_maxQ)
        self._set_data('detector_ADC_satQ', ADC_satQ)
        self._set_data('detector_ADC_bit', ADC_bit)
        self._set_data('detector_fpn_type', fpn_type)
        self._set_data('detector_fpn_offset', fpn_offset)
        self._set_data('detector_fpn_noise', fpn_noise)
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
        print '\tA/D Converter = ', self.detector_ADC_bit, 'bit'
        print '\tADC Max. Charge = ', self.detector_ADC_maxQ, 'electron'
        print '\tADC Sat. Charge = ', self.detector_ADC_satQ, 'electron'
        print '\tADC Fixed Pattern Noise :', self.detector_fpn_type
        print '\t\tOffset = ', self.detector_fpn_offset, 'ADC count'
        print '\t\tNoise  = ', self.detector_fpn_noise, 'ADC count'
        print '\tDark Current = ', self.detector_dark_current, 'electron/pixel/sec'
        print '\tReadout = ', self.detector_readout, 'electron'
        print '\tExcess	 = ', self.detector_excess
        print '\tEM gain = ', self.detector_emgain, 'x'
        print '\tQuantum Efficiency = ', self.detector_qeff[int(self.psf_wavelength)-int(self.wave_length[0])]



    def set_Movie(self, image_file_dir = None,
			movie_filename='./movies/movie.mp4', \
			cleanup_image_file_dir=False) :

	if image_file_dir is None:
	    image_file_dir = tempfile.mkdtemp(dir=os.getcwd())
	    cleanup_image_file_dir = True
            
	self._set_data('movie_image_file_dir', image_file_dir)
	self._set_data('movie_cleanup_image_file_dir', cleanup_image_file_dir)
	self._set_data('movie_filename', movie_filename)



    def set_Output(self, output_file_dir = None) :

        if output_file_dir is None:
            output_file_dir = tempfile.mkdtemp(dir=os.getcwd())

        self._set_data('output_file_dir', output_file_dir)



    def set_DataFile(self, hdf5_file_path_list, observable=None) :

	# read hdf5 lattice file
	for hdf5_file_path in hdf5_file_path_list :

	    try :

		hdf5_file = h5py.File(hdf5_file_path, 'r')
	
		species = hdf5_file['species']
		lattice = hdf5_file['lattice_info/HCP_group']
	        dataset = hdf5_file['data']
	
		# particle data in time-series
		data = []

                start = self.detector_start_time
                end   = self.detector_end_time

		for i in dataset :

		    data_i = dataset[i]
	            time = data_i.attrs['t']

		    if (time >= start and time < end) :

			particles = []

			for j in hdf5_file['data/'+str(i)+'/particles'] :
			    particles.append(j)

			element = [time, particles]
			data.append(element)

		data.sort(lambda x, y:cmp(x[0], y[0]))

		# get data
		self._set_data('spatiocyte_data', data)

		# get species properties
		self._set_data('spatiocyte_species_id', copy.copy(map(lambda x : x[0], species)))
		self._set_data('spatiocyte_index',      copy.copy(map(lambda x : x[1], species)))
		self._set_data('spatiocyte_diffusion',  copy.copy(map(lambda x : x[3], species)))
		self._set_data('spatiocyte_radius',     copy.copy(map(lambda x : x[2], species)))

		# get lattice properties
		self._set_data('spatiocyte_lattice_id', copy.copy(map(lambda x : x[0], lattice)))
		self._set_data('spatiocyte_lengths',    copy.copy(map(lambda x : x[1], lattice)))
		self._set_data('spatiocyte_VoxelRadius',   copy.copy(lattice[0][2]))
		self._set_data('spatiocyte_theNormalizedVoxelRadius', copy.copy(lattice[0][3]))
		self._set_data('spatiocyte_theStartCoord', copy.copy(lattice[0][4]))
		self._set_data('spatiocyte_theRowSize',    copy.copy(lattice[0][6]))
		self._set_data('spatiocyte_theLayerSize',  copy.copy(lattice[0][5]))
		self._set_data('spatiocyte_theColSize',    copy.copy(lattice[0][7]))

		hdf5_file.close()


            except Exception, e :
	                if not self.ignore_open_errors:
	                    raise
	                print 'Ignoring error: ', e

	# set observable
        if observable is None :
            index = [True for i in range(len(self.spatiocyte_index))]

        else :
            index = map(lambda x :  True if x.find(observable) > -1 else False, self.spatiocyte_index)

        self.spatiocyte_observables = copy.copy(index)



	# Visualization error	
	if self.spatiocyte_species_id is None:
	    raise VisualizerError('Cannot find species_id in any given hdf5 files')

	if len(self.spatiocyte_data) == 0:
	    raise VisualizerError('Cannot find particles dataset in any given hdf5 files: ' \
	                    + ', '.join(hdf5_file_path_list))

	if len(self.spatiocyte_index) == 0 :
	    raise VisualizerError('Cannot find lattice dataset in any given hdf5 files: ' \
	                    + ', '.join(hdf5_file_path_list))




    def set_BCFile(self, bc_file_path_list=None) :

	# read CSV Boundary condition file
	for bc_file_path in bc_file_path_list :

            try :

                csv_file = open(bc_file_path, 'r')

		dataset = []

		for row in csv.reader(csv_file) :
		    dataset.append(row)

                ### header
		header = dataset[0]

                ### data for boundary condition
                bc_data = []

		for j in range(1, len(dataset[1]), 3) :

		    c_id = [float(dataset[1][j]), float(dataset[1][j+1]), float(dataset[1][j+2])]
		    bc_data.append(c_id)
			
		bc_data.sort(lambda x, y:cmp(x[0], y[0]))

                # get data
                self._set_data('spatiocyte_bc_switch', True)
                self._set_data('spatiocyte_bc', bc_data)


            except Exception, e :
                        if not self.ignore_open_errors:
                            raise
                        print 'Ignoring error: ', e

    def set_efficiency(self, array, index=1):
        if len(array[0]) < 3:
            index = 1
        d = {int(item[0]): float(item[index]) for item in array}
        return [d.get(wl, 0.0) for wl in self.wave_length]

    def set_Optical_path(self) :

	# (1) Illumination path : Light source --> Sample
	self.set_Illumination_path()

	# (2) Detection path : Sample --> Detector
	self.set_Detection_path()



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
        f_s  = self.scanlens_focal_length
        f_t1 = self.tubelens_focal_length1

	w_tube = (f_t1/f_s)*w_BE

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

	self.source_flux = numpy.array(map(lambda x, y : (2*x)/(numpy.pi*y**2)*numpy.exp(-2*(r*1e-9/y)**2), N_z, w_z))


    def set_Detection_path(self) :

        wave_length = self.psf_wavelength*1e-9

	# Magnification
        Mag = 1.0

        # (2) objective lens
        f_obj = self.objective_focal_length

	# (3) tube lens
        f_t1 = self.tubelens_focal_length1

	# Magnification : Obj to tube lens
	Mag = (f_t1/f_obj)*Mag

        # (4) scan lens
        f_s = self.scanlens_focal_length

        # (5) focusing lens in front of detector
        f_t2 = self.tubelens_focal_length2

        # Magnification : scan to pinhole lens
        Mag = (f_t2/f_s)*Mag

	# set image scaling factor
        voxel_radius = self.spatiocyte_VoxelRadius

	view = self.detector_pixel_length/(2.0*voxel_radius)
	zoom = self.detector_zoom

	self.image_scaling = view/(Mag*zoom)
	self.image_magnification = Mag*zoom

	print 'Magnification :', Mag
	print 'Scaling :', self.image_scaling


        # Detector PSF
        self.set_PSF_detector()




    def set_PSF_detector(self) :

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

	# For normalization
	norm = map(lambda x : True if x > 1e-4 else False, I)


	# PSF : Fluorophore
	psf_fl = None

	if (self.fluorophore_type == 'Gaussian' or
	    self.fluorophore_type == 'Point-like' ) :

            I0 = 1.0
            Ir = sum(map(lambda x : x*numpy.exp(-0.5*(r/self.psf_width[0])**2), norm))
            Iz = sum(map(lambda x : x*numpy.exp(-0.5*(z/self.psf_width[1])**2), norm))

	    psf_fl = numpy.array(map(lambda x : I0*Ir*x, Iz))


	else :

	    # make the norm and wave_length array shorter
#	    psf_fl = 0
#
#	    for i in range(len(norm)) :
#
#		if norm[i] is True :
#	            psf_fl += self.get_PSF_fluorophore(r, z, wave_length[i])
#
#	    psf_fl = psf_fl/sum(norm)

	    psf_fl = numpy.sum(I)*self.get_PSF_fluorophore(r, z, wave_length)

	self.fluorophore_psf = psf_fl

    def get_PSF_fluorophore(self, r, z, wave_length):
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

        J0 = numpy.array(map(lambda x: j0(x*alpha*rho), r))
        Y = numpy.array(
            map(lambda x: 2*numpy.exp(-2*1.j*x*gamma*rho**2)*rho*drho, z))

        I = numpy.array(map(lambda x: x*J0, Y))
        I_sum = I.sum(axis=2)

        psf = numpy.array(map(lambda x: abs(x)**2, I_sum))
        Norm = numpy.amax(psf)

        return psf/Norm



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




class EPIFMVisualizer() :

	'''
	EPIFM Visualization class of e-cell simulator
	'''

	def __init__(self, configs=EPIFMConfigs(), effects=PhysicalEffects()) :

		assert isinstance(configs, EPIFMConfigs)
		self.configs = configs

                assert isinstance(effects, PhysicalEffects)
                self.effects = effects

		"""
		Check and create the folders for image and output files.
		"""
		if not os.path.exists(self.configs.movie_image_file_dir):
		    os.makedirs(self.configs.movie_image_file_dir)
		#else:
		#    for file in os.listdir(self.configs.movie_image_file_dir):
		#	os.remove(os.path.join(self.configs.movie_image_file_dir, file))

                if not os.path.exists(self.configs.output_file_dir):
                    os.makedirs(self.configs.output_file_dir)

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



#        def __del__(self):
#
#            if self.configs.movie_cleanup_image_file_dir :
#
#                for parent_dir, dirs, files in os.walk(self.configs.movie_image_file_dir, False) :
#                    for file in files :
#                        os.remove(os.path.join(parent_dir, file))
#
#                    os.rmdir(parent_dir)



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

        def polar2cartesian_coordinates(self, r, t, x, y):
            X, Y = numpy.meshgrid(x, y)

            new_r = numpy.sqrt(X*X + Y*Y)
            new_t = numpy.arctan2(X, Y)

            ir = interp1d(r, numpy.arange(len(r)), bounds_error=False)
            it = interp1d(t, numpy.arange(len(t)))

            new_ir = ir(new_r.ravel())
            new_it = it(new_t.ravel())

            new_ir[new_r.ravel() > r.max()] = len(r)-1
            new_ir[new_r.ravel() < r.min()] = 0

            return numpy.array([new_ir, new_it])

        def polar2cartesian(self, grid, coordinates, shape):
            r = shape[0] - 1
            psf_cart = numpy.empty([2 * r + 1, 2 * r + 1])
            psf_cart[r:, r:] = map_coordinates(grid, coordinates, order=0).reshape(shape)
            psf_cart[r:, :r] = psf_cart[r:, :r:-1]
            psf_cart[:r, :] = psf_cart[:r:-1, :]
            return psf_cart



	def get_depth_of_focus(self, p_i, p_0) :

		# get focal point
		x_0, y_0, z_0 = p_0

		# get particle position
		x_i, y_i, z_i = p_i

		wave_length = self.configs.psf_wavelength

		NA  = self.configs.objective_NA
		Mag = self.configs.image_magnification

		# constant
		a = self.effects.depth_of_focus_a
		b = self.effects.depth_of_focus_b

		# depth of focal : Berek's formula
		depth_1st = a/(NA**2)*wave_length
		depth_2nd = b/(Mag*NA)

		dof = depth_1st + depth_2nd

		# get particle depth from focal plane
		depth  = abs(x_i - x_0)

		if (depth < dof) :
		    depth = 0
		else :
		    depth = depth - dof

		return depth



	def get_intensity(self, time, pid, source_psf, source_max) :

		delta = self.configs.detector_exposure_time

                # linear conversion of intensity
                Ratio = self.effects.conversion_ratio
                intensity = Ratio*source_psf

		# Photobleaching process : Exponential decay
		if (self.effects.bleaching_switch == True) :

		    zeta  = self.effects.bleaching_rate
		    state = source_psf/source_max

		    self.effects.bleaching_state[pid] -= zeta*state*self.effects.bleaching_state[pid]*delta

		    #print pid, self.effects.bleaching_state[pid]

                    intensity = self.effects.bleaching_state[pid]*intensity


                # Slow-Blinking process : Power-law probability distribution
                if (self.effects.blinking_switch == True) :

		    state  = self.effects.blinking_state[pid]
		    period = self.effects.blinking_period[pid]

		    numpy.random.seed()
		    value = scipy.random.uniform(0, 1)
		    prob = self.effects.get_Prob_blinking(state, period)

		    if (value > prob) :
		        self.effects.blinking_state[pid]  = int(bool(self.effects.blinking_state[pid]-1))
		        self.effects.blinking_period[pid] = 0

		    else :
		        self.effects.blinking_period[pid] += delta

                    #print pid, value, prob, self.effects.blinking_state[pid], self.effects.blinking_period[pid]

                    intensity = self.effects.blinking_state[pid]*intensity


		return intensity

        def get_signal(self, time, pid, s_index, p_i, p_b, p_0, Norm):
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

            if (d_s < len(d)):
                source_depth = d_s
            else:
                source_depth = d[-1]

            # beam lateral position
            rr = numpy.sqrt((y_i-y_0)**2 + (z_i-z_0)**2)

            if (rr < len(r)):
                source_radius = rr
            else:
                source_radius = r[-1]

            # normalization
            unit_time = 1.0/self.configs.detector_fps

            radius = 1e-9

            unit_area = radius**2
            norm = unit_area*unit_time

            # get illumination PSF
            source_psf = norm*self.configs.source_flux[int(source_depth)][int(source_radius)]

            # signal conversion : Output Intensity = Physics * PSF (Beam)
            Ratio = self.effects.conversion_ratio
            Intensity = Ratio * source_psf

            # fluorophore axial position
            if self.effects.depth_of_focus_switch is True:
                d_f = self.get_depth_of_focus(p_i, p_0)
            else:
                d_f = abs(x_i - x_0)

            if (d_f < len(d)):
                fluo_depth = d_f
            else:
                fluo_depth = d[-1]

            # signal conversion : Output PSF = Intensity * PSF(Fluorophore)
            signal = Norm * Intensity * self.fluo_psf[int(fluo_depth)]

            return signal

        def get_noise(self, signal):
            # detector noise
            Nr = self.configs.detector_readout
            DC = self.configs.detector_dark_current
            T  = self.configs.detector_exposure_time
            Fn = numpy.sqrt(self.configs.detector_excess)
            M  = self.configs.detector_emgain

            # Noise calculation defined by HAMAMATSU Photonics
            sigma2 = Fn**2*M**2*(signal + DC*T)+ (Nr)**2
            noise  = numpy.sqrt(sigma2)

            return noise



        def set_boundary_plane(self, cell, p_b, p_0) :

                voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9

		# set focal point
		x_0, y_0, z_0 = p_0*voxel_size

                # set beam center
                x_b, y_b, z_b = p_b*voxel_size

		# data : boundary condition
		data = self.configs.spatiocyte_bc


                for i in range(len(data)) :

                    # particles relative coordinate
                    p_i = numpy.array(data[i])*voxel_size
		    x_i, y_i, z_i = p_i
		    #x_i, y_i, z_i = data[i]

		    # fluorophore axial position
		    #if (self.effects.depth_of_focus_switch == True) :
			#d_f = self.get_depth_of_focus(p_i, p_0)
		    #else :
			#d_f = abs(x_i - x_0)

		    #if (int(abs(x_i - x_0)) == 0) :
		    cell[int(z_i)][int(y_i)] = 1


		# Detector Output
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

		for i in range(Nw_cell/Np) :
		    for j in range(Nh_cell/Np) :

			# get signal
			signal = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])

			if (signal > 0) :
			    cell_pixel[i][j] = True
			else :
			    cell_pixel[i][j] = False


		return cell_pixel

        def overwrite_signal(self, cell, signal, p_i):
            # particle position
            x_i, y_i, z_i = p_i

            # z-axis
            Nz_cell = len(cell)
            Nz_signal = len(signal)
            Nr = len(self.configs.radial)

            z_to = z_i + Nr
            z_from = z_i - Nr

            if z_to > Nz_cell:
                dz_to = z_to - Nz_cell
                z0_to = int(Nz_cell)
                zi_to = int(Nz_signal - dz_to)
            else:
                dz_to = Nz_cell - (z_i + Nr)
                z0_to = int(Nz_cell - dz_to)
                zi_to = int(Nz_signal)

            if z_from < 0:
                dz_from = abs(z_from)
                z0_from = 0
                zi_from = int(dz_from)
            else:
                dz_from = z_from
                z0_from = int(dz_from)
                zi_from = 0

            ddz = (z0_to - z0_from) - (zi_to - zi_from)

            if ddz > 0:
                z0_to = z0_to - ddz
            if ddz < 0:
                zi_to = zi_to - ddz

            # y-axis
            Ny_cell = cell.size/Nz_cell
            Ny_signal = signal.size/Nz_signal

            y_to = y_i + Nr
            y_from = y_i - Nr

            if y_to > Ny_cell:
                dy_to = y_to - Ny_cell
                y0_to = int(Ny_cell)
                yi_to = int(Ny_signal - dy_to)
            else:
                dy_to = Ny_cell - (y_i + Nr)
                y0_to = int(Ny_cell - dy_to)
                yi_to = int(Ny_signal)

            if y_from < 0:
                dy_from = abs(y_from)
                y0_from = 0
                yi_from = int(dy_from)
            else:
                dy_from = y_from
                y0_from = int(dy_from)
                yi_from = 0

            ddy = (y0_to - y0_from) - (yi_to - yi_from)

            if ddy > 0:
                y0_to = y0_to - ddy
            if ddy < 0:
                yi_to = yi_to - ddy

            # add to cellular plane
            cell[z0_from:z0_to, y0_from:y0_to] += signal[zi_from:zi_to, yi_from:yi_to]

        def get_molecule_plane(self, cell, time, data, pid, p_b, p_0):
            voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9

            # particles coordinate, species and lattice IDs
            c_id, s_id, l_id = data

            sid_array = numpy.array(self.configs.spatiocyte_species_id)
            s_index = (numpy.abs(sid_array - int(s_id))).argmin()

            if self.configs.spatiocyte_observables[s_index] is True:
                # particles coordinate in real(nm) scale
                pos = self.get_coordinate(c_id)
                p_i = numpy.array(pos)*voxel_size

                # photons spreading to all directions
                Norm = 1.00/(4.0*numpy.pi)

                if self.configs.spatiocyte_index[s_index] == 'dR':
                    Norm = 2.00/(4.0*numpy.pi)

                # get signal matrix
                signal = self.get_signal(time, pid, s_index, p_i, p_b, p_0, Norm)

                # add signal matrix to image plane
                self.overwrite_signal(cell, signal, p_i)

        def output_frames(self, num_div=1):
            # set A/D Converter
            self.set_ADConverter()

            self.set_fluo_psf()

            start = self.configs.detector_start_time
            end = self.configs.detector_end_time
            exposure_time = self.configs.detector_exposure_time
            num_timesteps = int(math.ceil((end - start) / exposure_time))

            if self.use_multiprocess():
                num_processes = min(multiprocessing.cpu_count(), num_timesteps)
                n, m = divmod(num_timesteps, num_processes)
                # when 10 tasks is distributed to 4 processes,
                # number of tasks of each process must be [3, 3, 2, 2]
                chunks = [n + 1 if i < m else n for i in range(num_processes)]

                processes = []
                start_index = 0
                for chunk in chunks:
                    stop_index = start_index + chunk
                    process = multiprocessing.Process(
                        target=self.output_frames_each_process,
                        args=(start_index, stop_index))
                    processes.append(process)
                    start_index = stop_index
                for process in processes:
                    process.start()
                for process in processes:
                    process.join()
            else:
                self.output_frames_each_process(0, num_timesteps)

        def output_frames_each_process(self, start_count, stop_count):
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

            # set boundary condition
            if self.configs.spatiocyte_bc_switch is True:
                bc = numpy.zeros(shape=(Nz, Ny))
                bc = self.set_boundary_plane(bc, p_b, p_0)

            exposure_time = self.configs.detector_exposure_time

            detector_start_time = self.configs.detector_start_time
            time = detector_start_time + exposure_time * start_count

            # data-time interval
            t0 = self.configs.spatiocyte_data[0][0]
            t1 = self.configs.spatiocyte_data[1][0]

            delta_data = t1 - t0
            delta_time = int(round(exposure_time/delta_data))
            count0 = int(round(detector_start_time / exposure_time))

            for count in range(start_count + count0, stop_count + count0):
                # set image file name
                image_file_name = os.path.join(
                    self.configs.movie_image_file_dir,
                    self.configs.movie_image_file_name_format % (count))

                print 'time : ', time, ' sec (', count, ')'

                # define cell
                cell = numpy.zeros(shape=(Nz, Ny))

                count_start = (count - count0) * delta_time
                count_end = (count - count0 + 1) * delta_time

                frame_data = self.configs.spatiocyte_data[count_start:count_end]

                # loop for frame data
                for i, (i_time, data) in enumerate(frame_data):
                    print '\t', '%02d-th frame : ' % (i), i_time, ' sec'
                    # loop for particles
                    for j, data_j in enumerate(data):
                        self.get_molecule_plane(cell, i_time - time, data_j, j, p_0, p_b)

                if numpy.amax(cell) > 0:
                    if self.configs.spatiocyte_bc_switch:
                        camera = self.detector_output(cell, bc)
                    else:
                        camera = self.detector_output(cell)

                    camera.astype('uint%d' % (self.configs.detector_ADC_bit))

                    # save image to file
                    toimage(camera, cmin=1600, cmax=6000).save(image_file_name)
                time += exposure_time

	def overwrite_smeared(self, cell_pixel, photon_dist, i, j) :

		# i-th pixel
		Ni_pixel = len(cell_pixel)
		Ni_pe    = len(photon_dist)

                i_to   = i + Ni_pe/2
                i_from = i - Ni_pe/2

                if (i_to > Ni_pixel) :

                    di_to = i_to - Ni_pixel

                    i0_to = int(Ni_pixel)
                    i1_to = int(Ni_pe - di_to)

                else :

                    di_to = Ni_pixel - (i + Ni_pe/2)

                    i0_to = int(Ni_pixel - di_to)
                    i1_to = int(Ni_pe)

                if (i_from < 0) :

                    di_from = abs(i_from)

                    i0_from = 0
                    i1_from = int(di_from)

                else :

                    di_from = i_from

                    i0_from = int(di_from)
                    i1_from = 0

                ddi = (i0_to - i0_from) - (i1_to - i1_from)

                if (ddi > 0) : i0_to = i0_to - ddi
                if (ddi < 0) : i1_to = i1_to - ddi

                # j-th pixel
		Nj_pixel = len(cell_pixel[0])
		Nj_pe    = len(photon_dist[0])

                j_to   = j + Nj_pe/2
                j_from = j - Nj_pe/2

                if (j_to > Nj_pixel) :

                    dj_to = j_to - Nj_pixel

                    j0_to = int(Nj_pixel)
                    j1_to = int(Nj_pe - dj_to)

                else :

                    dj_to = Nj_pixel - (j + Nj_pe/2)

                    j0_to = int(Nj_pixel - dj_to)
                    j1_to = int(Nj_pe)

                if (j_from < 0) :

                    dj_from = abs(j_from)

                    j0_from = 0
                    j1_from = int(dj_from)

                else :

                    dj_from = j_from

                    j0_from = int(dj_from)
                    j1_from = 0

                ddj = (j0_to - j0_from) - (j1_to - j1_from)

                if (ddj > 0) : j0_to = j0_to - ddj
                if (ddj < 0) : j1_to = j1_to - ddj

		# add to cellular plane
                cell_pixel[i0_from:i0_to, j0_from:j0_to] += photon_dist[i1_from:i1_to, j1_from:j1_to]

		return cell_pixel

        def detector_output(self, cell, bc=None):
            # Detector Output
            voxel_radius = self.configs.spatiocyte_VoxelRadius
            voxel_size = (2.0*voxel_radius)/1e-9

            Nw_pixel = self.img_width
            Nh_pixel = self.img_height

            Np = int(self.configs.image_scaling*voxel_size)

            # image in nm-scale
            Nw_camera = Nw_pixel*Np
            Nh_camera = Nh_pixel*Np

            Nw_cell = len(cell)
            Nh_cell = len(cell[0])

            if Nw_camera > Nw_cell:
                w_cam_from = int((Nw_camera - Nw_cell)/2.0)
                w_cam_to = w_cam_from + Nw_cell
                w_cel_from = 0
                w_cel_to = Nw_cell
            else:
                w_cam_from = 0
                w_cam_to = Nw_camera
                w_cel_from = int((Nw_cell - Nw_camera)/2.0)
                w_cel_to = w_cel_from + Nw_camera

            if Nh_camera > Nh_cell:
                h_cam_from = int((Nh_camera - Nh_cell)/2.0)
                h_cam_to = h_cam_from + Nh_cell
                h_cel_from = 0
                h_cel_to = Nh_cell
            else:
                h_cam_from = 0
                h_cam_to = int(Nh_camera)
                h_cel_from = int((Nh_cell - Nh_camera)/2.0)
                h_cel_to = h_cel_from + Nh_camera

            # image in nm-scale
            plane = cell[w_cel_from:w_cel_to, h_cel_from:h_cel_to]

            # convert image in nm-scale to pixel-scale
            cell_pixel = numpy.zeros(shape=(Nw_cell/Np, Nh_cell/Np))

            # Photon distribution
            for i in range(Nw_cell/Np):
                for j in range(Nh_cell/Np):

                    # get photons
                    photons = numpy.sum(plane[i*Np:(i+1)*Np,j*Np:(j+1)*Np])

                    if photons > 0:

                        # get crosstalk
                        if self.effects.detector_crosstalk_switch is True:

                            width = self.effects.detector_crosstalk_width

                            n_i = numpy.random.normal(0, width, photons)
                            n_j = numpy.random.normal(0, width, photons)

                            #i_bins = int(numpy.amax(n_i) - numpy.amin(n_i))
                            #j_bins = int(numpy.amax(n_j) - numpy.amin(n_j))

                            smeared_photons, edge_i, edge_j = numpy.histogram2d(n_i, n_j, bins=(24, 24), range=[[-12,12],[-12,12]])

                            # smeared photon distributions
                            cell_pixel = self.overwrite_smeared(cell_pixel, smeared_photons, i, j)

                        else:

                            cell_pixel[i][j] = photons

            index = int(self.configs.psf_wavelength) - int(self.configs.wave_length[0])
            QE = self.configs.detector_qeff[index]

            # Photoelectron and ADC count distribution
            for i in range(Nw_cell/Np):
                for j in range(Nh_cell/Np):

                    # pixel position
                    pixel = (i, j)

                    # Detector : Quantum Efficiency
                    index = int(self.configs.psf_wavelength) - int(self.configs.wave_length[0])
                    QE = self.configs.detector_qeff[index]

                    # get signal (photoelectrons)
                    signal = QE*cell_pixel[i][j]

                    # get constant background (photoelectrons)
                    if self.effects.background_switch is True:
                        mean = self.effects.background_mean
                        background = QE*numpy.random.poisson(mean, None)
                    else:
                        background = 0

                    # get detector noise (photoelectrons)
                    noise = self.get_noise(signal + background)

                    # get EM Gain
                    M  = self.configs.detector_emgain

                    # get signal + background (photoelectrons)
                    PE = numpy.random.normal(M*(signal + background), noise, None)
                    #PE = M*(signal + background)

                    # A/D converter : Photoelectrons --> ADC counts
                    ADC = self.get_ADC_value(pixel, PE)
                    #ADC = PE

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
                    noise = self.get_noise(M*(signal + background))

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



	def set_ADConverter(self) :

	    Nw_pixel = self.img_width
	    Nh_pixel = self.img_height

	    # set ADC offset
	    ADC0  = self.configs.detector_fpn_offset
	    noise = self.configs.detector_fpn_noise

	    # set Fixed-Pattern noise
            FPN_type = self.configs.detector_fpn_type

	    if (FPN_type == None) :
                offset = numpy.array([ADC0 for i in range(Nw_pixel*Nh_pixel)])

	    if  (FPN_type == 'pixel') :
		offset = numpy.rint(numpy.random.normal(ADC0, noise, Nw_pixel*Nh_pixel))

	    elif (FPN_type == 'column') :
                column = numpy.random.normal(ADC0, noise, Nh_pixel)
                temporal = numpy.array([column for i in range(Nw_pixel)])

		offset = numpy.rint(temporal.reshape(Nh_pixel*Nw_pixel))

	    # set ADC gain
            ADC_bit  = self.configs.detector_ADC_bit
            ADC_satQ = self.configs.detector_ADC_satQ

            gain = numpy.array(map(lambda x : (ADC_satQ - 0)/(2**ADC_bit - x), offset))

	    print 'A/D Converter :', ADC_bit, '-bit'
	    print '\toffset =', ADC0
	    print '\tgain  =', (ADC_satQ - 0)/(2**ADC_bit - ADC0)
	    print '\t', '%s-FPN : sigma =' % (FPN_type), noise, 'ADC count'

            self.configs.detector_ADC_offset = offset.reshape([Nw_pixel, Nh_pixel])
            self.configs.detector_ADC_gain = gain.reshape([Nw_pixel, Nh_pixel])

        def get_ADC_value(self, pixel, photo_electron):
            # pixel position
            i, j = pixel

            # check non-linearity
            Q_max = self.configs.detector_ADC_satQ

            if photo_electron > Q_max:
                photo_electron = Q_max

            # convert photoelectron to ADC counts (Grayscale)
            k = self.configs.detector_ADC_gain[i][j]
            ADC0 = self.configs.detector_ADC_offset[i][j]
            ADC_max = 2**self.configs.detector_ADC_bit - 1

            ADC = photo_electron/k + ADC0

            if ADC > ADC_max:
                ADC = ADC_max

            if ADC < 0:
                ADC = 0

            return int(ADC)

        def make_movie(self) :
            """
            Make a movie by FFmpeg
            Requirement : Install FFmpeg (http://ffmpeg.org/)
            """
            input_image_filename = os.path.join(self.configs.movie_image_file_dir, self.configs.movie_image_file_name_format)
    
            # Set FFMPEG command
            ffmpeg_command  = 'ffmpeg -sameq -r "%s"' % (str(int(self.configs.detector_fps)))
            ffmpeg_command += ' -y -i "./%s/%s" ' % (self.configs.movie_image_file_dir, self.configs.movie_image_file_name_format)
            ffmpeg_command += self.configs.movie_filename
            #ffmpeg_command += ('" -vcodec rawvideo -pix_fmt yuyv422 ' + self._movie_filename)
    
            os.system(ffmpeg_command)



	def output_movie(self, num_div=1):
	    """
	    Output movie
	    """
	    self.output_frames(num_div=num_div)
	    self.make_movie()

        def use_multiprocess(self):
            envname = 'ECELL_MICROSCOPE_SINGLE_PROCESS'
            return not (envname in os.environ and os.environ[envname])

        def set_fluo_psf(self):
            depths = set()

            voxel_size = 2.0*self.configs.spatiocyte_VoxelRadius/1e-9
            start = self.configs.detector_start_time
            end = self.configs.detector_end_time
            exposure_time = self.configs.detector_exposure_time
            num_timesteps = int(math.ceil((end - start) / exposure_time))

            t0 = self.configs.spatiocyte_data[0][0]
            t1 = self.configs.spatiocyte_data[1][0]
            delta_data = t1 - t0
            delta_time = int(round(exposure_time / delta_data))

            for count in range(num_timesteps):
                Nz = int(self.configs.spatiocyte_lengths[0][2] * voxel_size)
                Ny = int(self.configs.spatiocyte_lengths[0][1] * voxel_size)
                Nx = int(self.configs.spatiocyte_lengths[0][0] * voxel_size)

                # focal point
                p_0 = numpy.array([Nx, Ny, Nz])*self.configs.detector_focal_point

                count_start = count * delta_time
                count_end = (count + 1) * delta_time

                frame_data = self.configs.spatiocyte_data[count_start:count_end]
                for _, data in frame_data:
                    for data_j in data:
                        c_id = data_j[0]
                        pos = self.get_coordinate(c_id)
                        p_i = numpy.array(pos) * voxel_size

                        # fluorophore axial position
                        d = self.configs.depth
                        if self.effects.depth_of_focus_switch is True:
                            d_f = self.get_depth_of_focus(p_i, p_0)
                        else:
                            d_f = abs(p_i[0] - p_0[0])
                        if d_f < len(d):
                            fluo_depth = d_f
                        else:
                            fluo_depth = d[-1]
                        depths.add(int(fluo_depth))

            depths = list(depths)
            if self.use_multiprocess():
                self.fluo_psf = {}
                num_processes = min(multiprocessing.cpu_count(), len(depths))
                n, m = divmod(len(depths), num_processes)
                chunks = [n + 1 if i < m else n for i in range(num_processes)]

                processes = []
                start_index = 0
                for chunk in chunks:
                    stop_index = start_index + chunk
                    p, c = multiprocessing.Pipe()
                    process = multiprocessing.Process(
                        target=self.get_fluo_psf_process,
                        args=(depths[start_index:stop_index], c))
                    processes.append((p, process))
                    start_index = stop_index
                for _, process in processes:
                    process.start()
                for pipe, process in processes:
                    self.fluo_psf.update(pipe.recv())
                    process.join()
            else:
                self.fluo_psf = self.get_fluo_psf(depths)

        def get_fluo_psf(self, depths):
            r = self.configs.radial
            theta = numpy.linspace(0, 90, 91)
            z = numpy.linspace(0, +r[-1], len(r))
            y = numpy.linspace(0, +r[-1], len(r))
            coordinates = self.polar2cartesian_coordinates(r, theta, z, y)

            psf_t = numpy.ones_like(theta)

            result = {}
            for depth in depths:
                psf_r = self.configs.fluorophore_psf[depth]

                psf_polar = numpy.array(map(lambda x: psf_t * x, psf_r))

                result[depth] = self.polar2cartesian(psf_polar, coordinates, (len(r), len(r)))
            return result

        def get_fluo_psf_process(self, depths, pipe):
            pipe.send(self.get_fluo_psf(depths))
