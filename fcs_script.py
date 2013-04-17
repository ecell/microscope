"""
    fcs_script.py:

    User script to create the image from the simulated Fluoroscence Correlation Spectroscopy (FCS)
"""

from pFcs_handler import FCSConfigs, FCSVisualizer

def test_fcs() :

	# create TIRF Microscopy
	fcs = FCSConfigs()

	fcs.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=500e-6, radius=0.32e-3)
	fcs.set_BeamExpander(expander_type='Keplerian', focal_length1=30e-3, focal_length2=20e-3, pinhole_radius=13e-6)
	fcs.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	#fcs.set_Fluorophore(fluorophore_type='Gaussian', wave_length=578, width=(70.0, 140.0))
	#fcs.set_Fluorophore(fluorophore_type='Point-like', wave_length=600)
	fcs.set_Objective(NA=1.20, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
	#fcs.set_DichroicMirror('FF562-Di03-25x36')
	#fcs.set_EmissionFilter('FF01-593_40-25')
	fcs.set_TubeLens1(focal_length=160e-3)
	fcs.set_ScanLens(focal_length=50e-3)
	fcs.set_PinholeLens(focal_length=200e-3, radius=15e-6)
	fcs.set_Detector(detector='EMCCD', zoom=10, emgain=200, focal_point=(0.5,0.5,0.5), \
			start_time=0, end_time=1, fps=1e+3, exposure_time=1e-3)
	fcs.set_Movie(image_file_dir='./images_fcs', movie_filename='./movies/fcs_movie.mp4')
	fcs.set_DataFile(['./data/lattice/test_cube1.h5'])

	# create image and movie
	create = FCSVisualizer(configs=fcs, output_file='./test_fcs.dat')
	#create.get_plots(plot_filename='./plots/fcs_plots.pdf')
	#create.output_frames(num_div=1)
	create.output_movie(num_div=1)


if __name__ == "__main__":

	test_fcs()

