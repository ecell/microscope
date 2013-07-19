"""
    egfm_script.py:

    User script to create the image from the simulated Confocal Microscopy (Wide-field view FCS)
"""
import sys
import os

from pFcs_handler import FCSConfigs, FCSVisualizer
from tirfm_handler  import TIRFMConfigs, TIRFMVisualizer
from effects_handler  import PhysicalEffects

def test_confm(t0, t1, Molecule) :

	# create Confocal Microscopy
	confm = FCSConfigs()

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
	confm.set_Movie(image_file_dir='./images_egfm_%s' % (Molecule), movie_filename='./movies/egfm_movie.mp4')
	confm.set_DataFile(['/home/kiwamoto/share/EGF_model_20130603/LogHDF5_ERK_1.h5'], observable=Molecule)

	# create image and movie
	create = FCSVisualizer(configs=confm)
	#create.get_plots(plot_filename='./plots/epifm.plots.pdf')
	create.output_frames(num_div=16)
	#create.output_movie(num_div=16)


def test_tirfm(t0, t1, Molecule) :

        # create TIRF Microscopy
        tirfm = TIRFMConfigs()

        tirfm.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=1e-3, radius=0.32e-3)
        #tirfm.set_BeamExpander(expander_type='Keplerian', focal_length1=300e-3, focal_length2=20e-3, pinhole_radius=23e-6)
        #tirfm.set_Fluorophore(fluorophore_type='mCherry')
        tirfm.set_Fluorophore(fluorophore_type='EGFP')
        #tirfm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=578, width=(70.0, 140.0))
        #tirfm.set_Fluorophore(fluorophore_type='Point-like', wave_length=578)
        tirfm.set_Objective(NA=1.49, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
        tirfm.set_EvanescentField(depth=200)
        tirfm.set_DichroicMirror('FF562-Di03-25x36')
        tirfm.set_EmissionFilter('FF01-593_40-25')
        tirfm.set_TubeLens1(focal_length=160e-3)
        tirfm.set_ScanLens(focal_length=50e-3)
        tirfm.set_TubeLens2(focal_length=200e-3)
        tirfm.set_Detector(detector='EMCCD', zoom=1, emgain=300, focal_point=(0.0,0.5,0.5), \
                        start_time=t0, end_time=t1, fps=1/33e-3, exposure_time=33e-3)
        tirfm.set_Movie(image_file_dir='./images_egfm_%s' % (Molecule), movie_filename='./movies/egfm_movie.mp4')
        tirfm.set_DataFile(['/home/kiwamoto/share/EGF_model_20130603/LogHDF5_Raf_1.h5'], observable=Molecule)
        tirfm.set_BCFile(['/home/kiwamoto/share/EGF_model_20130603/LogCoordinate_Shape.csv'])

        # create physical effects
        physics = PhysicalEffects()
        physics.set_Conversion(ratio=1e-6)
        #physics.set_Background(mean=30)
        #physics.set_DepthOfFocus(a=0.5, b=0.0)
        #physics.set_DetectorCrosstalk(width=2.0)
        #physics.set_Photobleaching(rate=1/5.0)
        #physics.set_Photoblinking(P0=0.5, alpha_on=0.0048, alpha_off=0.0055)

        # create image and movie
        create = TIRFMVisualizer(configs=tirfm, effects=physics)
        #create.get_plots(plot_filename='./plots/tirfm.plots.pdf')
        create.output_frames(num_div=16)
        #create.output_movie(num_div=16)



if __name__ == "__main__":

        t0 = float(sys.argv[1])
        t1 = float(sys.argv[2])

	#test_confm(t0, t1, "ERKPP")
	#test_confm(t0, t1, "ERK")
	test_tirfm(t0, t1, "Raf")

