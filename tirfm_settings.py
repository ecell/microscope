"""
 tirf_settigs.py:

    Settings for Fluoresence Microscope

	Fluorophore
	  -- Selected Fluorophore
	  -- Gaussian Point Spreading Function (PSF)
	Incident Beam Condition
	Excitation Filter
	Objective
        Dichroic Mirror
	Emission Filter
	Tube Lens
	Camera
	Movie

"""

import numpy
from rgb_colors import *

#-----------------------------
# General 
#-----------------------------
ignore_open_errors = False

#-----------------------------
# Particles list
#-----------------------------
particles_index = None
particles_list  = None

#-----------------------------
# Fluorophore
#-----------------------------
fluorophore_type = 'Gauss'
fluorophore_wavelength = numpy.array([i for i in range(300, 1500)])
fluorophore_wavenumber = numpy.array([2.*numpy.pi/fluorophore_wavelength[i] for i in range(len(fluorophore_wavelength))])
fluoex_eff  = numpy.array([0.0 for i in range(len(fluorophore_wavelength))])
fluoem_eff  = numpy.array([0.0 for i in range(len(fluorophore_wavelength))])

fluorophore_radial = numpy.array([1.0*i for i in range(1000)])
fluorophore_depth  = numpy.array([1.0*i for i in range(1000)])
fluorophore_signal = numpy.array([[0.0 for i in range(len(fluorophore_radial))] for j in range(len(fluorophore_depth))])
fluorophore_rgb    = numpy.array([(0, 0, 0) for j in range(len(fluorophore_depth))])

#-----------------------------
# PSF
#-----------------------------
psf_wavelength = 600 # nm
psf_width  = (200, 200)	# Gaussian function (radial width, lateral width) [nm]
psf_cutoff = (400, 100)	# cutoff range (radius, depth)
psf_file_name_format = 'psf_%04d.png'	# Image file name

penetration_depth = 200 # nm

#-----------------------------
# Incident Beam Condition
#----------------------------- 
beam_switch  = False       
beam_wavelength = 600.  # nm
beam_intensity = 3.108e-3 # W/cm^2
beam_pinhole_radius = len(fluorophore_radial)

#-----------------------------
# Mirror
#-----------------------------
mirror_position   = 0.8 # epi when the postion is at 0

#-----------------------------
# Excitation Filter
#-----------------------------
excitation_switch = False
excitation_eff = numpy.array([0.0 for i in range(len(fluorophore_wavelength))])

#-----------------------------
# Objective
#-----------------------------
objective_switch = False
objective_NA  = 1.45
objective_Ng  = 1.52
objective_Nm  = 1.37
objective_mag = 60
objective_efficiency   = 0.90
objective_sin_alpha    = None
objective_sin_critical = None
objective_eff = numpy.array([0.0 for i in range(len(fluorophore_wavelength))])
objective_glassthickness = 3.0 # mm

#-----------------------------
# Dichroic Mirror
#-----------------------------
dichroic_switch = False
dichroic_eff = numpy.array([0.0 for i in range(len(fluorophore_wavelength))])

#-----------------------------
# Emission Filter
#-----------------------------
emission_switch = False
emission_eff = numpy.array([0.0 for i in range(len(fluorophore_wavelength))])

#-----------------------------
# Tube Lens
#-----------------------------
tubelens_switch = False
tubelens_mag = 4
tubelens_range = (380, 780)  # nm
tubelens_efficiency = 0.9987

#-----------------------------
# Pinhole
#-----------------------------
pinhole_radius = len(fluorophore_radial)

#-----------------------------
# Camera
#-----------------------------
camera_switch = False
camera_type = 'Perfect'
camera_base_position = (-2.0, 0.5, 0.5) # Base position of x,y,z (This unit is world_size)
camera_focal_point   = ( 0.0, 0.5, 0.5) # Focal point of x,y,z (This unit is world_size)
camera_image_size   = (512, 512)        # camera image size in pixels
camera_pixel_length = 16.0e-6           # Pixel size in micro-m scale
camera_zoom = 1.0                       # Zoom Zoom-in > 1.0 > Zoom-out
camera_gain = 1.0                       # Zoom Zoom-in > 1.0 > Zoom-out
camera_start_time = 0.0
camera_end_time = None
camera_exposure_time  = 0.033
camera_fps = 1 #5-30
camera_sat_charge = 30000.
camera_ADC_bit   = 16
camera_ADC_const = 1.00
camera_ADC_offset = 2000

camera_readout = 0.0
camera_dark_current = 0.0
camera_excess = 1.0
camera_emgain = 1.0

camera_blue  = numpy.array([1.0 for i in range(len(fluorophore_wavelength))])
camera_green = numpy.array([1.0 for i in range(len(fluorophore_wavelength))])
camera_red   = numpy.array([1.0 for i in range(len(fluorophore_wavelength))])

#-----------------------------
# Image/Movie
#-----------------------------
movie_background_color = RGB_LIGHT_SLATE_GRAY
movie_image_file_dir = "./images"
movie_image_file_name_format = 'image_%04d.png' # Must be compatible with FFmpeg's input-file notation
movie_cleanup_image_file_dir = False
movie_filename = "./movies/movie.mp4"

#-----------------------------
# Spatiocyte
#-----------------------------
spatiocyte_data = []

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
# SNR
#-----------------------------
photon_number = numpy.array([i+1 for i in range(10000)])
absolute_snr  = numpy.array([[0.0 for i in range(len(photon_number))] for j in range(len(fluorophore_wavelength))])
relative_snr  = numpy.array([[1.0 for i in range(len(photon_number))] for j in range(len(fluorophore_wavelength))])
ideal_snr = numpy.array([photon_number[i]/numpy.sqrt(photon_number[i]) for i in range(len(photon_number))])
ideal_relsnr  = numpy.array([1.0 for i in range(len(photon_number))])


