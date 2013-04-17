"""
    kinesin_script.py:

    User script to create the image from the simulated Epifluoroscence Microscopy (only for off-lattice model)
"""

from kinesin_handler import KinesinConfigs, KinesinVisualizer

def test_kinesin() :

	# create EPIF Microscopy for off-lattice model (Kinesin)
	epifm = KinesinConfigs()

	epifm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=10e-3, radius=0.32e-3)
	epifm.set_BeamExpander(expander_type='Keplerian', focal_length1=300e-3, focal_length2=20e-3, pinhole_radius=23e-6)
	epifm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	#epifm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=578, width=(70.0, 140.0))
	#epifm.set_Fluorophore(fluorophore_type='Point-like', wave_length=578)
	epifm.set_Objective(NA=1.49, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
	epifm.set_DichroicMirror('FF562-Di03-25x36')
	epifm.set_EmissionFilter('FF01-593_40-25')
	epifm.set_TubeLens1(focal_length=160e-3)
	epifm.set_ScanLens(focal_length=50e-3)
	epifm.set_TubeLens2(focal_length=200e-3)
	epifm.set_Detector(detector='EMCCD', zoom=2, emgain=1200, focal_point=(0.5,0.0,0.0), \
			start_time=0, end_time=16, fps=1/30e-3, exposure_time=30e-3)
	epifm.set_Movie(image_file_dir='./images_kinesin', movie_filename='./movies/kinesin_movie.mp4')
	epifm.set_DataFile(['./data/lattice/satyaKinesinVisualLog.csv'])

	# create image and movie
	create = KinesinVisualizer(configs=epifm)
	#create.get_plots(plot_filename='./plots/kinesin_plots.pdf')
	#create.output_frames(num_div=16)
	create.output_movie(num_div=16)

	############################################################################
	#   Note 
	#	you may need to rotate focal position
	#	change two things;
	#		(1) get_coordinate() ---- exchange point_x <--> point_y
	#		(2) output_frames() ----- exchange Nx <--> Ny
        ############################################################################



if __name__ == "__main__":

	test_kinesin()

