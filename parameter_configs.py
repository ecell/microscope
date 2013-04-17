"""
 param_configs.py:

    Parameter Configurations for Fluoresence Microscope

	Fluorophore
	  -- Selected Fluorophore (Airy)
	  -- Gaussian Point Spreading Function (PSF)
	Light Source
          -- Selected Light Source
          -- Gaussian Beam
	Mirror Position
	Excitation Filter
	Objective
        Dichroic Mirror
	Emission Filter
	Scan Lens
	Tube Lens
	Detector

"""

import numpy

#-----------------------------
# General 
#-----------------------------
ignore_open_errors = False

radial = numpy.array([1.0*i for i in range(1000)])
depth  = numpy.array([1.0*i for i in range(1000)])
wave_length = numpy.array([i for i in range(300, 1000)])
wave_number = numpy.array([2.*numpy.pi/wave_length[i] for i in range(len(wave_length))])

image_scaling = 1.0

#-----------------------------
# Particles list
#-----------------------------
particles_index = None
particles_list  = None

#-----------------------------
# Fluorophore
#-----------------------------
fluorophore_type = 'Gauss'
fluorophore_lifetime = 1e-4 # sec
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

penetration_depth = 200 # nm

#-----------------------------
# Scattering Matrix of photon-matter interactions 
#-----------------------------
scattering_matrix = None

#-----------------------------
# Incident Beam Condition
#----------------------------- 
source_switch  = False       
source_type = 'LASER'
source_wavemode = 'TEM00'
source_M2factor = 1.0
source_wavelength = 600. # nm
source_power = 20e-3 # W
source_radius = 0.32e-3 # m
source_intensity  = 3.108e-3 # W/cm^2
source_divergence = numpy.array([0.0 for i in range(len(depth))])
source_flux = numpy.array([[0.0 for i in range(len(radial))] for j in range(len(depth))])

#-----------------------------
# Beam Expander
#-----------------------------
expander_type = 'Keplerian' # Keplerian or Gallilean
expander_pinhole_radius = 23e-3 # m
expander_focal_length1 = 5e-3  # m
expander_focal_length2 = 25e-3 # m

#-----------------------------
# Mirror
#-----------------------------
mirror_position = 0.0

#-----------------------------
# Excitation Filter
#-----------------------------
excitation_switch = False
excitation_eff = numpy.array([0.0 for i in range(len(wave_length))])

#-----------------------------
# Objective
#-----------------------------
objective_switch = False
objective_NA  = 1.45
objective_Ng  = 1.52
objective_Nm  = 1.37
objective_focal_length = 1.9e-3 # m
objective_efficiency   = 0.90
objective_sin_alpha    = None
objective_sin_critical = None
objective_eff = numpy.array([0.0 for i in range(len(wave_length))])
objective_glassthickness = 3.0e-3 # m

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
# Tube Lens
#-----------------------------
tubelens_switch = False
tubelens_focal_length1 = 160e-3 # m
tubelens_focal_length2 = 200e-3 # m

#-----------------------------
# Scan Lens
#-----------------------------
tubelens_switch = False
tubelens_focal_length = 50e-3 # m

#-----------------------------
# Pinhole
#-----------------------------
pinholelens_switch = False
pinholelens_focal_length = 50e-3 # m
pinholelens_radius = 23e-6 # m

#-----------------------------
# Detector (CCD camera, PMT, ADP, .... etc)
#-----------------------------
detector_switch = False
detector_type = 'Perfect'
detector_base_position = (-2.0, 0.5, 0.5) # Base position of x,y,z (This unit is world_size)
detector_focal_point   = ( 0.0, 0.5, 0.5) # Focal point of x,y,z (This unit is world_size)
detector_image_size   = (512, 512)        # detector image size in pixels
detector_pixel_length = 16.0e-6           # Pixel size in micro-m scale
detector_zoom = 1.0                       # Zoom Zoom-in > 1.0 > Zoom-out
detector_gain = 1.0                       # Zoom Zoom-in > 1.0 > Zoom-out
detector_start_time = 0.0
detector_end_time = None
detector_exposure_time  = 0.033
detector_fps = 1 #5-30
detector_sat_charge = 370000.
detector_max_charge = 600000.
detector_ADC_bit   = 16
detector_ADC_const = 5.8
detector_ADC_offset = 2000

detector_readout = 0.0
detector_dark_current = 0.0
detector_excess = 1.0
detector_emgain = 1.0
detector_background = 0.0

detector_qeff  = numpy.array([1.0 for i in range(len(wave_length))])
#detector_blue  = numpy.array([1.0 for i in range(len(wave_length))])
#detector_green = numpy.array([1.0 for i in range(len(wave_length))])
#detector_red   = numpy.array([1.0 for i in range(len(wave_length))])

#-----------------------------
# Image/Movie
#-----------------------------
movie_background_color = (0, 0, 0)
movie_image_file_dir = "./images"
movie_image_file_name_format = 'image_%04d.png' # Must be compatible with FFmpeg's input-file notation
#movie_image_file_name_format = 'image_%04d.tiff' # Must be compatible with FFmpeg's input-file notation
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

