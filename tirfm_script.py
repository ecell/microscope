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
        tirfm.set_BeamExpander(expander_type='Keplerian', focal_length1=100e-3, focal_length2=20e-3, pinhole_radius=23e-6)
	tirfm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	#tirfm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=650, width=(70.0, 140.0))
	#tirfm.set_Fluorophore(fluorophore_type='Point-like', wave_length=600)
	tirfm.set_Mirror(position=0.9)
	tirfm.set_Objective(NA=1.49, Nm=1.37, efficiency=0.90)
	#tirfm.set_DichroicMirror('FF562-Di03-25x36')
	#tirfm.set_EmissionFilter('FF01-593_40-25')
        tirfm.set_TubeLens(focal_length=160e-3)
        tirfm.set_ScanLens(focal_length=50e-3)
	tirfm.set_Detector(detector='EMCCD', zoom=20.0, focal_point=(0.0,0.5,0.5), \
			start_time=0, end_time=1, fps=1e+3, exposure_time=1e-3)
	tirfm.set_Movie(image_file_dir='./images_tirfm', movie_filename='./movies/tirfm_movie.mp4')
	#tirfm.set_DataFile(['./data/lattice/test_model.h5'])
	tirfm.set_DataFile(['./data/lattice/test_cube1.h5'])

	# create image and movie
	create = TIRFMVisualizer(configs=tirfm)
	#create.get_plots(plot_filename='./plots/tirfm_plots.pdf', wave_length=610)
	create.output_frames(num_div=16)
	#create.output_movie(num_div=16)



if __name__ == "__main__":

	test_tirfm()

