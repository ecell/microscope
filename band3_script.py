"""
    band3_script.py:

    User script to create the image from the simulated Epifluoroscence Microscopy (EPIFM)
"""
import sys
import os

#from epifm_handler import EPIFMConfigs, EPIFMVisualizer
from kinesin_handler import KinesinConfigs, KinesinVisualizer

def test_b3c(t0, t1) :

	# create EPIF Microscopy
	epifm = KinesinConfigs()

	epifm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=10e-3, radius=0.32e-3)
	epifm.set_BeamExpander(expander_type='Keplerian', focal_length1=300e-3, focal_length2=20e-3, pinhole_radius=23e-6)
	#epifm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	epifm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=578, width=(10.0, 20.0))
	#epifm.set_Fluorophore(fluorophore_type='Point-like', wave_length=578)
	epifm.set_Objective(NA=1.49, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
	#epifm.set_DichroicMirror('FF562-Di03-25x36')
	#epifm.set_EmissionFilter('FF01-593_40-25')
	epifm.set_TubeLens1(focal_length=160e-3)
	epifm.set_ScanLens(focal_length=50e-3)
	epifm.set_TubeLens2(focal_length=200e-3)
	epifm.set_Detector(detector='EMCCD', zoom=1, emgain=100, pixel_length=0.16e-6, focal_point=(0.0,0.5,0.5), \
			start_time=t0, end_time=t1, fps=1, exposure_time=1)
	epifm.set_Movie(image_file_dir='./images_b3c', movie_filename='./movies/band3_cluster.mp4')
	epifm.set_DataFile(['./data/lattice/band3_cluster.csv'])

	# create image and movie
	create = KinesinVisualizer(configs=epifm)
	#create.get_plots(plot_filename='./plots/epifm_plots.pdf')
	create.output_frames(num_div=16)
	#create.output_movie(num_div=16)


if __name__ == "__main__":

	t0 = float(sys.argv[1])
	t1 = float(sys.argv[2])

	test_b3c(t0, t1)


