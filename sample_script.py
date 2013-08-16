"""
    tirfm_script.py:

    User script to create the image from the simulated
    Total Internal Refractive Fluoroscence Microscopy (TIRFM)
"""
import sys
from tirfm_handler import TIRFMConfigs, TIRFMVisualizer
from effects_handler import PhysicalEffects


def test_tirfm(t0, t1):
    # create TIRF Microscopy
    tirfm = TIRFMConfigs()
    tirfm.set_LightSource(source_type='LASER', wave_mode='TEM00',
                          M2_factor=1.0, wave_length=473,
                          power=1.0e-3, radius=10e-6)
    tirfm.set_Fluorophore(fluorophore_type='Tetramethylrhodamine(TRITC)')
    tirfm.set_EvanescentField(depth=200)
    tirfm.set_Objective(NA=1.49, Nm=1.37,
                        focal_length=1.90e-3, efficiency=0.90)
    tirfm.set_TubeLens1(focal_length=160e-3)
    tirfm.set_ScanLens(focal_length=50e-3)
    tirfm.set_TubeLens2(focal_length=200e-3)
    tirfm.set_Detector(detector='EMCCD', zoom=1, emgain=300,
                       focal_point=(0.0, 0.5, 0.5),
                       start_time=t0, end_time=t1,
                       fps=1.0/3.3e-3, exposure_time=33e-3)
    tirfm.set_Movie(image_file_dir='./images',
                    movie_filename='./movies/tirfm_movie.mp4')
    tirfm.set_DataFile(['./data/lattice/test_model_11.h5'])

    # create physical effects
    physics = PhysicalEffects()
    physics.set_Conversion(ratio=1e-6)
    physics.set_DetectorCrosstalk(width=2.00)

    # create image and movie
    create = TIRFMVisualizer(configs=tirfm, effects=physics)
    create.output_frames(num_div=16)


if __name__ == "__main__":
    t0 = float(sys.argv[1])
    t1 = float(sys.argv[2])

    test_tirfm(t0, t1)
