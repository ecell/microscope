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
        
from scipy.special import j0
from scipy.integrate import quad

def integrand(rho, r, z, wave_length) :

	NA = 1.49
        k = 2.0*scipy.pi/wave_length
        alpha = k*NA
        gamma = k*(NA/2)**2

	func = 2*j0(r*alpha*rho)*scipy.exp(-2*1.j*z*gamma*rho**2)*rho

	return func


def psfint(r, z, wave_length) :

	return abs(quad(integrand, 0, 1, args=(r, z, wave_length))[0])**2


def get_PSF(r, z, wave_length) :
    fluorophore_psf = numpy.zeros((1000,600))

    NA = 1.49#self.objective_NA
    N = 50#100
    drho = 1.0/N
    rho = numpy.linspace(1.0/N, 1, 50)

    k = 2.0*numpy.pi/wave_length
    alpha = k*NA
    gamma = k*(NA/2)**2

    J0 = numpy.array(map(lambda y : map(lambda x : j0(x*y*rho), r), alpha))
    Y  = numpy.array(map(lambda y : map(lambda x : 2*numpy.exp(-2*1.j*x*y*rho**2)*rho*drho, z), gamma))

    for i in range(len(wave_length)) :
        I  = numpy.array(map(lambda x : x*J0[i], Y[i]))
        I_sum = I.sum(axis=2)
        I_abs = map(lambda x : abs(x)**2, I_sum)

        fluorophore_psf += I_abs

        print wave_length[i], I_abs[0][0]


if __name__ == "__main__":

	r = numpy.linspace(0,599,600)
	z = numpy.linspace(0,999,1000)
	wave_length = numpy.linspace(300,999,700)

	result = map(lambda x : psfint(x, 0, 600), r)
	print result
	#get_PSF(r, z, wave_length)

