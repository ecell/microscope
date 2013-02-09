"""
    tirfm_script.py:

    User script to create the image from the simulated TIRF Microscopy
"""

from tirfm_handler import TIRFMSettings, TIRFMVisualizer

def test_fcs() :

	# create TIRF Microscopy
	fcs = TIRFMSettings()

	# Spectral Arrangement
	fcs.set_IncidentBeam(wave_length=473, intensity=50.00, excitation=None)
	#fcs.set_Fluorophore(fluorophore_type='Qdot 605')
	fcs.set_Fluorophore(fluorophore_type='Gaussian', wave_length=650, width=(70.0, 140.0))
	#fcs.set_Fluorophore(fluorophore_type='Point-like', wave_length=600)
	fcs.set_Mirror(position=0.0)
	fcs.set_Objective(mag=60, NA=1.49, Nm=1.37, efficiency=0.90)
	#fcs.set_DichroicMirror('FF562-Di03-25x36')
	#fcs.set_EmissionFilter('FF01-593_40-25')
	fcs.set_TubeLens(mag=4)
	fcs.set_Pinhole(radius=87)
	fcs.set_Detector(detector='EMCCD', zoom=1.0, focal_point=(0.5,0.5,0.5), \
			start_time=0, end_time=1, fps=1e+3, exposure_time=1e-3)
	#fcs.set_Movie(image_file_dir='./images', movie_filename='./movies/test_model.mp4')
	fcs.set_DataFile(['./data/lattice/test_cube0.h5'])

	# create image and movie
	create = TIRFMVisualizer(settings=fcs, output_file='./data/fcs/test_cube0.dat')
	#create.get_plots(plot_filename='./plots/spectral_plots.pdf', wave_length=610)
	create.output_frames(num_div=16)
	#create.output_movie(num_div=16)



if __name__ == "__main__":

	test_fcs()

