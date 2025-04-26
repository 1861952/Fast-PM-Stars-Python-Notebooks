# General imports:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Astropy imports:
from astropy.io import fits
from astropy import units as u
from astropy import constants as const

#Scipy imports
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
from scipy.stats import powerlaw

# Astropy imports 
from astropy.modeling.models import BlackBody
from astropy.visualization import quantity_support

# qp import
import qp

pi = np.pi
# Function for getting radial distribution
def get_radial_distribution(sample_size=100000, max_radius = 2000):
    # raised to the power of a-1, so squared
    a = 3
    mean, var, skew, kurt = powerlaw.stats(a, moments='mvsk')    
    # distance-distribution
    radial_dist = powerlaw(a, scale=max_radius)
    return radial_dist

# Function for getting distribution using qp model
def get_qp_distribution(bin_edges, counts, sample_size = 100000):
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    midpoints = 0.5*(bin_edges[:-1]+bin_edges[1:])
    # already normalized when check_input is True
    pdfs = counts
    
    # qp.Ensemble class has objects representing a set of probability density functions
    ens_h = qp.Ensemble(qp.hist, data=dict(bins=bin_edges, pdfs=pdfs))
    return ens_h

# Function for making f_lambda spline using Phoenix model
# f_lambda function using table
def f_lambda_func(wavelength_array):
    # Use 1e8 instead of 10**8 for non-integers
    return (f_lambda*1e8)*(wavelength_array*1e-8)**2/(const.c.cgs.value)
    
    

# F_nu function(with f_nu as function of wavelength, transmission, over normal integral)
def passband_flux(flux, transmissionArr, wavelengthsArr):
    int_n = integrate.trapezoid(flux*transmissionArr/wavelengthsArr.cgs.value, x=wavelengthsArr.cgs.value)
    int_d = integrate.trapezoid(transmissionArr/wavelengthsArr.cgs.value, x=wavelengthsArr.cgs.value)
    if int_d==0:
        print(int_d)
    return int_n/int_d



# Returns and compares calculated values of magnitude
# f_nu (surface flux), R (radius of star), d (distance from observer to star), and F_0 (zero-point flux) against given "actual" magnitude

def get_and_compare_mag(f_nu, R, d, F_0, actual_mag=None):
    if actual_mag is None:
        F_M = (f_nu * 4*pi*R**2)/(4*pi*d**2)
        # Absolute magnitude
        M = -2.5 * np.log10(F_M/F_0)
        return M
        

    else:
        F_M = (f_nu * 4*pi*R**2)/(4*pi*d**2)
        # Absolute magnitude
        M = -2.5 * np.log10(F_M/F_0)
        print(f"Calculated: {M}, Actual: {actual_mag}")
        return M