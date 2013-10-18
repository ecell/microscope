"""
    pointscan_confm_script.py:

    User script to create the image from the simulated Point-scanning Confocal Microscopy
"""
import sys
import os

from pointscan_confm_handler import PointScanConfocalConfigs, PointScanConfocalVisualizer
from effects_handler import PhysicalEffects

def test_confm(t0, t1) :

	# create Point-scanning Confocal Microscopy
	confm = PointScanConfocalConfigs()

        confm.set_LightSource(source_type='LASER', wave_length=532, power=100e-6, radius=200e-9)
	confm.set_Fluorophore(fluorophore_type='EGFP')
	#confm.set_Fluorophore(fluorophore_type='Gaussian', wave_length=508, width=(70.0, 140.0))
	#confm.set_DichroicMirror('FF562-Di03-25x36')
	#confm.set_EmissionFilter('FF01-593_40-25')
	confm.set_Pinhole(radius=8e-6)
	confm.set_Magnification(Mag=336)
	confm.set_Detector(detector='PMT', zoom=1, emgain=1e+6, focal_point=(0.4,0.5,0.5), \
			exposure_time=0.50, bandwidth=50e+3, mode='Pulse')
	confm.set_ADConverter(bit=16, offset=2000, fullwell=60000)
	confm.set_OutputData(image_file_dir='./images_ERK_3')
	confm.set_InputData('/home/kiwamoto/ecell3/ecell3-spatiocyte/Data/EGF_model_10/scale3_1/LogCoordinate_ERK_3', start=t0, end=t1)

	# create physical effects
	physics = PhysicalEffects()
	physics.set_Conversion(ratio=1e-6)
	#physics.set_Background(mean=30)
	#physics.set_DetectorCrosstalk(width=1.00)

	# create image and movie
	create = PointScanConfocalVisualizer(configs=confm, effects=physics)
	create.output_frames(num_div=16)



if __name__ == "__main__":

        t0 = float(sys.argv[1])
        t1 = float(sys.argv[2])

	test_confm(t0, t1)

