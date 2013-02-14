import sys
import os
import copy
import tempfile
import time 
import math
import operator
import random
        
import scipy
import numpy 
import pylab
        

def do_loop() :

	r = numpy.array([1.00*i for i in range(600)])

	z = numpy.linspace(-r[-1], +r[-1], 2*len(r)-1)
        y = numpy.linspace(-r[-1], +r[-1], 2*len(r)-1)

        Z, Y = numpy.meshgrid(z, y)
        R = numpy.sqrt(Z**2 + Y**2)

        # get PSF (Beam)
        beam_psf = copy.copy(R)
        beam_array = numpy.array([numpy.exp(-0.5*(r[i]/70.)**2) for i in range(len(r))])

        for i in range (len(R)) :
            for j in range(len(R[i])) :

                rr = R[i][j]

                if (rr < len(beam_array)) :
                    beam_psf[i][j] = beam_array[int(rr)]
                else :
                    beam_psf[i][j] = beam_array[-1]


	fig = pylab.figure()
	spec_scale = numpy.linspace(numpy.amin(beam_psf), numpy.amax(beam_psf), 200, endpoint=True)
	pylab.contour(Z, Y,  beam_psf, spec_scale, linewidth=0.1, color='k')
	pylab.contourf(Z, Y, beam_psf, cmap=pylab.cm.jet)
	pylab.show()




if __name__ == "__main__":

	do_loop()


