"""
    egfm_script.py:

    User script to create the image from the simulated Confocal Microscopy (Wide-field view FCS)
"""
import sys
import os

from confm_handler import ConfocalConfigs, ConfocalVisualizer

def test_confm(t0, t1) :

	# create Confocal Microscopy
	confm = ConfocalConfigs()

        confm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=500e-6, radius=0.32e-3)
        confm.set_BeamExpander(expander_type='Keplerian', focal_length1=30e-3, focal_length2=20e-3, pinhole_radius=13e-6)
	confm.set_Fluorophore(fluorophore_type='EGFP')
	#confm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=508, width=(70.0, 140.0))
	#confm.set_Fluorophore(fluorophore_type='Point-like', wave_length=508)
        confm.set_Objective(NA=1.20, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
        confm.set_TubeLens1(focal_length=160e-3)
        confm.set_ScanLens(focal_length=50e-3)
        confm.set_PinholeLens(focal_length=200e-3, radius=13e-6)
	confm.set_Detector(detector='EMCCD', zoom=1, emgain=5, focal_point=(0.5,0.5,0.5), \
			start_time=t0, end_time=t1, fps=1, exposure_time=1)
	confm.set_Movie(image_file_dir='./images_egfm_ERK', movie_filename='./movies/egfm_movie.mp4')
	confm.set_DataFile(['/home/kiwamoto/share/EGF_model_20130423/vlog_EGFmodel_15.h5'], observable="ERK")

	# create image and movie
	create = ConfocalVisualizer(configs=confm)
	#create.get_plots(plot_filename='./plots/epifm.plots.pdf')
	create.output_frames(num_div=16)
	#create.output_movie(num_div=16)



if __name__ == "__main__":

        t0 = float(sys.argv[1])
        t1 = float(sys.argv[2])

	test_confm(t0, t1)

