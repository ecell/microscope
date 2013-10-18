"""
    tirfm_script.py:

    User script to create the image from the simulated TIRF Microscope (TIRFM)
"""
import sys
import os

from tirfm_handler   import TIRFMConfigs, TIRFMVisualizer
from effects_handler import PhysicalEffects

def test_tirfm(t0, t1) :

	# create TIRF Microscopy
	tirfm = TIRFMConfigs()
	tirfm.set_LightSource(source_type='LASER', wave_length=473, power=1.0e-3, radius=10e-6)
	tirfm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
	#tirfm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=578, width=(20.0, 20.0))
	tirfm.set_DichroicMirror('FF562-Di03-25x36')
	tirfm.set_EmissionFilter('FF01-593_40-25')
	tirfm.set_Magnification(Mag=336)
	tirfm.set_Detector(detector='EMCCD', zoom=1, emgain=900, focal_point=(0.0,0.5,0.5), exposure_time=30e-3)
	tirfm.set_ADConverter(bit=16, offset=2000, fullwell=370000) # for EMCCD
	#tirfm.set_ADConverter(bit=16, offset=2000, fullwell=30000) # for CMOS
	#tirfm.set_OutputData(image_file_dir='./images_dicty_01_tirfm_x100')
	tirfm.set_OutputData(image_file_dir='./images_test')
	#tirfm.set_InputData('/home/masaki/ecell3/latest/data/csv/simple_dicty_00', start=t0, end=t1)
	tirfm.set_InputData('/home/masaki/ecell3/latest/data/csv/test', start=t0, end=t1)

	# create physical effects
	physics = PhysicalEffects()
	physics.set_Conversion(ratio=1e-6)
	#physics.set_Background(mean=30)
	#physics.set_DetectorCrosstalk(width=1.00)
	#physics.set_Photobleaching(rate=1/2.00)
	#physics.set_Photoblinking(P0=0.5, alpha_on=0.0048, alpha_off=0.0055)

	# create image and movie
	create = TIRFMVisualizer(configs=tirfm, effects=physics)
	create.output_frames(num_div=16)



if __name__ == "__main__":

        t0 = float(sys.argv[1])
        t1 = float(sys.argv[2])

	test_tirfm(t0, t1)


