import numpy as np
import pandas as pd

from astropy.io import fits
from astropy import units as u
from astropy import constants as const

#Scipy imports
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import integrate
from scipy.stats import powerlaw

from pathlib import Path


# Reads in Ascii file and returns wavelengths in cm and transmissions
def ascii_to_array(filename):
    # Read the text file into a DataFrame
    data = pd.read_csv(filename, sep=' ', header=None, names=['Wavelength', 'Transmission'])

    # wavelength in Angstrom
    column1 = (data['Wavelength'].values)*u.AA
    # Probability of photon going through
    column2 = data['Transmission'].values
    return column1, column2


# Reads in Fits file and returns spline for flux at surface (Phoenix Model)

def fits_to_spline(fileURL, log_g):
    hdul = fits.open(fileURL)
    data_fits = hdul[1].data
    # Wavelengths from fits table
    wavLen = data_fits['WAVELENGTH']
    # f_lamda values for specific log(g) 
    f_lambda = data_fits[log_g]
    # Makes spline based table
    # convert from Angstrom to cm
    wavelengthArr_cm = 1e-8*wavLen
    f_lambda_table_per_cm = f_lambda*1e8  
    # Returns spline at surface based on Phoenix model
    # convert from f_lambda to f_nu
    return InterpolatedUnivariateSpline(wavelengthArr_cm, f_lambda_table_per_cm*wavelengthArr_cm**2/const.c.cgs.value, k=1, ext="zeros")


# Finds nearest value (absolute value) in an array
def find_nearest_value(target, array):
    target = float(target)
    diff_arr = np.abs(np.asarray(array)-target)
    index = (np.abs(array - target)).argmin()
    return array[index]
# Creates a dictionary where key is (T_eff, log(g)) and value is the file name

def create_dict(foldername, glob_str):
    # Get folder
    path = Path(foldername)
    path.glob(glob_str)

    # create dictionary : key is (T, log(g)) value is 'filename'
    file_dict = dict()
    for filepath in path.glob("sp_*"):
        filename = filepath.name
        fr = filepath.name.removeprefix("sp_t").removesuffix("nc_m-0.5.txt")
        T_eff = float(fr[0:4])
        g_mks = float(fr[5:])
        log_g_cgs = np.log10(g_mks*float(100))
        file_dict.update({(T_eff, log_g_cgs): filename})
        
    # np.unique automatically sorts
    T_eff_arr = np.unique([temp for temp, _lgg in file_dict])
    log_g_arr = np.unique([lgg for _temp, lgg in file_dict])
    T_eff_arr = np.asarray(T_eff_arr)
    log_g_arr = np.asarray(log_g_arr)
    return T_eff_arr, log_g_arr, file_dict


# Finds correct Sonora model file name based on T_eff and log(g) input
def find_file(T_eff_actual, log_g_actual, foldername, glob_str):
    T_eff_array, log_g_array, file_dict = create_dict(foldername, glob_str) 
    closest_T_eff = find_nearest_value(T_eff_actual, T_eff_array)
    closest_log_g = find_nearest_value(log_g_actual, log_g_array)
    g_mks = np.round(float(10)**(closest_log_g-float(2)))
    filename = f"{foldername}/sp_t{int(closest_T_eff)}g{int(g_mks)}nc_m-0.5.txt" 
    return filename

# Read in Sonora Model Spectra Files
def sonora_to_spline(filename):
    data = pd.read_csv(filename, skiprows=3, sep='\s+', names=['Wavelength', 'Flux'])
    # UNABLE TO DO ASTROPY UNITS FOR BOTH
    # convert from microns to centimeters
    wavelengths = data['Wavelength']*1e-4
    wavelengths = wavelengths[::-1]
    # observer flux
    Fluxes = data['Flux']
    Fluxes = Fluxes[::-1]
    # Makes spline based on table
    # Returns spline at surface based on Phoenix model
    # (already in f_nu so don't need to convert)
    return InterpolatedUnivariateSpline(wavelengths, Fluxes, k=1, ext="zeros")
