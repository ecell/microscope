"""
    kinesin_script.py:

    User script to create the image from the simulated Epifluoroscence Microscopy (only for off-lattice model)
"""

from kinesin_handler import KinesinConfigs, KinesinVisualizer

def test_kinesin() :

	# create Confocal Microscopy for off-lattice model (Kinesin)
	confm = KinesinConfigs()

        confm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=500e-6, radius=0.32e-3)
        confm.set_BeamExpander(expander_type='Keplerian', focal_length1=30e-3, focal_length2=20e-3, pinhole_radius=13e-6)
        confm.set_Fluorophore(fluorophore_type='EGFP')
        #confm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=508, width=(70.0, 140.0))
        #confm.set_Fluorophore(fluorophore_type='Point-like', wave_length=508)
        confm.set_Objective(NA=1.20, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
        confm.set_TubeLens1(focal_length=160e-3)
        confm.set_ScanLens(focal_length=50e-3)
        confm.set_PinholeLens(focal_length=200e-3, radius=13e-6)
	confm.set_Detector(detector='PMT', zoom=2, emgain=1e+7, focal_point=(0.5,0.5,0.5), \
			start_time=0, end_time=30, fps=1/30e-3, exposure_time=30e-3)
	confm.set_Movie(image_file_dir='./images_kinesin', movie_filename='./movies/kinesin_m100_movie.mp4')
	#confm.set_DataFile(['./data/lattice/satyaKinesinVisualLog.csv'])
	confm.set_DataFile(['./data/lattice/CoordinateLog_act0.15_m100.csv'])

	# create image and movie
	create = KinesinVisualizer(configs=confm)
	#create.get_plots(plot_filename='./plots/kinesin_plots.pdf')
	#create.output_frames(num_div=16)
	create.output_movie(num_div=16)

	############################################################################
	#   Note
	#
	#	you may need to rotate focal position
	#	change following things;
	#		(1) get_coordinate() ---- exchange point_x <--> point_y
	#		(2) output_frames() ----- exchange Nx <--> Ny
	#
        ############################################################################



if __name__ == "__main__":

	test_kinesin()

