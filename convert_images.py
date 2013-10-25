import sys
import math
import copy

import numpy
from scipy.misc import imread, toimage

cmin = 2000
cmax = 2010
#cmax = 12000

def convert(file_in, file_out) :

	i = 0

	max_count = 0

	while (True) :
	    try :
		input_image  = imread(file_in + '/image_%07d.png' % (i))
	    except Exception :
		break

	    output_image = file_out + '/image_%07d.png' % (i)

	    amax = numpy.amax(input_image)
	    amin = numpy.amin(input_image)

	    if (max_count < amax) :
		max_count = amax

	    print i, amax, amin

	    toimage(input_image, cmin=cmin, cmax=cmax).save(output_image)

	    i += 1

	print 'Max count : ', max_count, 'ADC'


if __name__=='__main__':

	#file_in  = '/home/masaki/microscopy/images_dicty_02_epifm_zaxis03'
	#file_in  = '/home/masaki/microscopy/images_dicty_02_line_zaxis09'
	#file_in  = '/home/masaki/microscopy/images_test'
	file_in  = '/home/masaki/microscopy/images_dicty_01_point'
	file_out = '/home/masaki/microscopy/images'

	convert(file_in, file_out)

