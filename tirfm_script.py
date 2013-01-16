"""
    tirfm_script.py:

    User script to create the image from the simulated TIRF Microscopy
"""

from tirfm_handler import TIRFMSettings, TIRFMVisualizer

def test_image() :

	# create TIRF Microscopy
	tirfm = TIRFMSettings()

	# Particles
	#tirfm.set_Particles(['Raf', 'Raf_RasGTP', 'RafA', 'RafA_Pase1', 'MEK_RafA', 'MEKP_RafA'])

	# Spectral Arrangement
	tirfm.set_IncidentBeam(wave_length=590, intensity=1e+2, excitation=None)
	#tirfm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	tirfm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=650, width=(70.0, 140.0))
	#tirfm.set_Fluorophore(fluorophore_type='Point-like', wave_length=600)
	tirfm.set_Mirror(position=0.9)
	tirfm.set_Objective(mag=60, NA=1.49, Nm=1.37, efficiency=0.90)
	#tirfm.set_DichroicMirror('FF562-Di03-25x36')
	#tirfm.set_EmissionFilter('FF01-593_40-25')
	tirfm.set_TubeLens(mag=4)
	#tirfm.set_Pinhole(radius=150)
	tirfm.set_Camera(camera='EMCCD', zoom=3.0, focal_point=(0.0,0.5,0.5), \
			start_time=0, end_time=33, fps=30, exposure_time=0.033)
	tirfm.set_Movie(image_file_dir='./images', movie_filename='./movies/test_model.mp4')
	tirfm.set_DataFile(['./data/lattice/test_model_volume.h5'])
	#tirfm.set_DataFile(['./data/lattice/test_cube.h5'])

	# create image and movie
	create = TIRFMVisualizer(settings=tirfm)
	#create.get_plots(plot_filename='./plots/spectral_plots.pdf', wave_length=610)
	create.output_movie(num_div=16)



if __name__ == "__main__":

	test_image()

