"""
    param_configs.py:

    Parameter Configurations for Fluoresence Microscope

	Fluorophore (PSF)
	  -- Selected Fluorophore (Airy)
	  -- Gaussian Point Spreading Function
	Light Source
	Excitation Filter
        Dichroic Mirror
	Emission Filter
	Pinhole/Slit
	Detector

"""

import numpy

#-----------------------------
# General 
#-----------------------------
ignore_open_errors = False
electron_charge = 1.602e-19 # C

radial = numpy.array([1.0*i for i in range(1000)])
depth  = numpy.array([1.0*i for i in range(1000)])
wave_length = numpy.array([i for i in range(300, 1000)])
wave_number = numpy.array([2.*numpy.pi/wave_length[i] for i in range(len(wave_length))])

#-----------------------------
# Fluorophore
#-----------------------------
fluorophore_type = 'Gauss'
fluorophore_lifetime = 10e-9 # sec
#fluoex_eff  = numpy.array([0.0 for i in range(len(wave_length))])
#fluoem_eff  = numpy.array([0.0 for i in range(len(wave_length))])
fluoex_eff  = [0.0 for i in range(len(wave_length))]
fluoem_eff  = [0.0 for i in range(len(wave_length))]

fluorophore_psf = numpy.array([[0.0 for i in range(len(radial))] for j in range(len(depth))])
#fluorophore_rgb = numpy.array([(0, 0, 0) for j in range(len(depth))])

#-----------------------------
# Fluorophore PSF 
#-----------------------------
psf_wavelength = 600 # nm
psf_width  = (200, 200)	# Gaussian function (radial width, lateral width) [nm]
psf_cutoff = (400, 100)	# cutoff range (radius, depth)
psf_file_name_format = 'psf_%04d.png'	# Image file name

#-----------------------------
# Illumination
#----------------------------- 
source_switch  = False       
source_type = 'LASER'
source_wavelength = 600. # nm
source_power = 20e-3 # W
source_radius = 20e-6 # m
source_depth  = 20e-6 # m
source_flux = numpy.array([[0.0 for i in range(len(radial))] for j in range(len(depth))])

#-----------------------------
# Resolution/Magnification
#-----------------------------
image_resolution = 16e-8
image_magnification = 100

#-----------------------------
# Excitation Filter
#-----------------------------
excitation_switch = False
excitation_eff = numpy.array([0.0 for i in range(len(wave_length))])

#-----------------------------
# Dichroic Mirror
#-----------------------------
dichroic_switch = False
dichroic_eff = numpy.array([0.0 for i in range(len(wave_length))])

#-----------------------------
# Emission Filter
#-----------------------------
emission_switch = False
emission_eff = numpy.array([0.0 for i in range(len(wave_length))])

#-----------------------------
# Pinhole
#-----------------------------
pinhole_radius = 16e-6 # m

#-----------------------------
# Slit
#-----------------------------
slit_size = 37e-6 # m

#-----------------------------
# Detector
#-----------------------------
detector_switch = False
detector_type = 'Perfect'
detector_base_position = (-2.0, 0.5, 0.5) # Base position of x,y,z (This unit is world_size)
detector_focal_point   = ( 0.0, 0.5, 0.5) # Focal point of x,y,z (This unit is world_size)
detector_image_size   = (512, 512)        # detector image size in pixels
detector_pixel_length = 16.0e-6           # Pixel size in micro-m scale
detector_zoom = 1.0                       # Zoom Zoom-in > 1.0 > Zoom-out
detector_gain = 1.0                       # Zoom Zoom-in > 1.0 > Zoom-out
detector_exposure_time  = 0.033
detector_bandwidth = 0
detector_mode = 'Pulse'
detector_readout = 0.0
detector_dark_current = 1.0e-12
detector_excess = 1.0
detector_emgain = 1.0
detector_background = 0.0

detector_qeff  = numpy.array([1.0 for i in range(len(wave_length))])
#detector_blue  = numpy.array([1.0 for i in range(len(wave_length))])
#detector_green = numpy.array([1.0 for i in range(len(wave_length))])
#detector_red   = numpy.array([1.0 for i in range(len(wave_length))])

#-----------------------------
# A/D Converter
#-----------------------------
ADConverter_fullwell = 370000.
ADConverter_bit  = 16
ADConverter_gain = 5.8
ADConverter_offset = 2000
ADConverter_fpn_type = None
ADConverter_fpn_count = 0.0

#-----------------------------
# Image
#-----------------------------
image_file_dir = "./images"
image_file_name_format = 'image_%07d.png'
image_file_cleanup_dir = False

#-----------------------------
# Spatiocyte
#-----------------------------
spatiocyte_start_time = 0
spatiocyte_start_end  = 1
spatiocyte_interval = 1e-3
spatiocyte_data = []
spatiocyte_observable = []

spatiocyte_species_id = []
spatiocyte_index      = []
spatiocyte_diffusion  = []
spatiocyte_radius     = []
spatiocyte_lattice_id = []
spatiocyte_lengths    = []

spatiocyte_VoxelRadius = 10e-9
spatiocyte_theNormalizedVoxelRadius = 0.5
spatiocyte_theStartCoord = None
spatiocyte_theRowSize    = None
spatiocyte_theLayerSize  = None
spatiocyte_theColSize    = None

#-----------------------------
# Spatiocyte Boundary condition
#-----------------------------
spatiocyte_bc_switch = False
spatiocyte_bc = []
