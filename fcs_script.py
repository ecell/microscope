"""
    fcs_script.py:

    User script to create the image from the simulated Fluoroscence Correlation Spectroscopy (FCS)
"""
import sys
#from fcs_handler import FCSConfigs, FCSVisualizer
from pFcs_handler import FCSConfigs, FCSVisualizer

def test_fcs(t0, t1, psf) :

	# create TIRF Microscopy
	fcs = FCSConfigs()

	fcs.set_LightSource(source_type='LASER', wave_mode='TEM00', M2_factor=1.00, wave_length=473, power=500e-6, radius=0.32e-3)
	fcs.set_BeamExpander(expander_type='Keplerian', focal_length1=30e-3, focal_length2=20e-3, pinhole_radius=13e-6)

	if (psf == 'airy') :
	    fcs.set_Fluorophore(fluorophore_type='Qdot 605')
	elif (psf == 'gauss') :
	    fcs.set_Fluorophore(fluorophore_type='Gaussian', wave_length=605, width=(100.0, 200.0))

	fcs.set_Objective(NA=1.20, Nm=1.37, focal_length=1.9e-3, efficiency=0.90)
	#fcs.set_DichroicMirror('FF562-Di03-25x36')
	#fcs.set_EmissionFilter('FF01-593_40-25')
	fcs.set_TubeLens1(focal_length=160e-3)
	fcs.set_ScanLens(focal_length=50e-3)
	fcs.set_PinholeLens(focal_length=200e-3, radius=13e-6)
	fcs.set_Detector(detector='PMT', zoom=2, emgain=1e+7, focal_point=(0.5,0.5,0.5), \
			start_time=t0, end_time=t1, fps=1/10e-5, exposure_time=10e-5)
	fcs.set_Movie(image_file_dir='./images_fcs_0100_%s' % (psf), movie_filename='./movies/fcs_movie.mp4')
	fcs.set_Output(output_file_dir='./output_fcs_0100_%s' % (psf))
	fcs.set_DataFile(['./data/lattice/test_fcs_0100.h5'])

	# create image and movie
	create = FCSVisualizer(configs=fcs)
	#create.get_plots(plot_filename='./plots/fcs_airy.pdf')
	create.output_frames(num_div=1)
	#create.output_movie(num_div=1)


if __name__ == "__main__":

	t0  = float(sys.argv[1])
	t1  = float(sys.argv[2])
	psf = str(sys.argv[3])

	test_fcs(t0, t1, psf)

