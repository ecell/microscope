"""
    tirfm_script.py:

    User script to create the image from the simulated Total Internal Refractive Fluoroscence  Microscopy (TIRFM)
"""

from tirfm_handler import TIRFMConfigs, TIRFMVisualizer

def test_tirfm() :

	# create TIRF Microscopy
	tirfm = TIRFMConfigs()

	# Spectral Arrangement
        tirfm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.0, wave_length=473, power=10e-3, radius=0.32e-3)
        tirfm.set_BeamExpander(expander_type='Keplerian', focal_length1=300e-3, focal_length2=20e-3, pinhole_radius=23e-6)
	tirfm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	#tirfm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=578, width=(70.0, 140.0))
	#tirfm.set_Fluorophore(fluorophore_type='Point-like', wave_length=600)
	tirfm.set_Mirror(position=0.9013)
	tirfm.set_Objective(NA=1.49, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
	tirfm.set_DichroicMirror('FF562-Di03-25x36')
	tirfm.set_EmissionFilter('FF01-593_40-25')
        tirfm.set_TubeLens1(focal_length=160e-3)
        tirfm.set_ScanLens(focal_length=50e-3)
        tirfm.set_TubeLens2(focal_length=200e-3)
	tirfm.set_Detector(detector='EMCCD', zoom=2, emgain=100, focal_point=(0.0,0.5,0.5), \
			start_time=10, end_time=10.033, fps=30, exposure_time=33e-3)
	#tirfm.set_Movie(image_file_dir='./images_tirfm', movie_filename='./movies/tirfm_movie.mp4')
	tirfm.set_DataFile(['./data/lattice/test_cube2.h5'])

	# create image and movie
	create = TIRFMVisualizer(configs=tirfm)
	create.get_plots(plot_filename='./plots/tirfm_plots.pdf')
	#create.output_frames(num_div=16)
	#create.output_movie(num_div=16)



if __name__ == "__main__":

	test_tirfm()

