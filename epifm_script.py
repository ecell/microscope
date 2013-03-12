"""
    epifm_script.py:

    User script to create the image from the simulated Epifluoroscence Microscopy (EPIFM)
"""

from epifm_handler import EPIFMConfigs, EPIFMVisualizer

def test_epifm() :

	# create EPIF Microscopy
	epifm = EPIFMConfigs()

	epifm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.0, wave_length=473, power=500e-6, radius=0.32e-3)
	#epifm.set_BeamExpander(expander_type='Keplerian', focal_length1=30e-3, focal_length2=40e-3, pinhole_radius=23e-6)
	epifm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	epifm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=473, width=(70.0, 140.0))
	#epifm.set_Fluorophore(fluorophore_type='Point-like', wave_length=600)
	epifm.set_Objective(NA=1.49, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
	#epifm.set_DichroicMirror('FF562-Di03-25x36')
	#epifm.set_EmissionFilter('FF01-593_40-25')
	epifm.set_TubeLens(focal_length=120e-3)
	epifm.set_ScanLens(focal_length=50e-3)
	epifm.set_Detector(detector='EMCCD', zoom=20.0, focal_point=(0.5,0.5,0.5), \
			start_time=0, end_time=1, fps=1e+4, exposure_time=1e-4)
	epifm.set_Movie(image_file_dir='./images', movie_filename='./movies/test_epifm.mp4')
	#epifm.set_DataFile(['./data/lattice/test_model.h5'])
	epifm.set_DataFile(['./data/lattice/test_cube2.h5'])

	# create image and movie
	create = EPIFMVisualizer(configs=epifm)
	#create.get_plots(plot_filename='./plots/spectral_plots.pdf', wave_length=610)
	create.output_frames(num_div=16)
	#create.output_movie(num_div=16)



if __name__ == "__main__":

	test_epifm()

