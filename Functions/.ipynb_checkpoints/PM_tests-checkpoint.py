import numpy as np
import pandas as pd
from scipy.stats import norm, powerlaw
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import Galactocentric, galactocentric_frame_defaults, ICRS, SkyCoord, Galactic
from Functions.PM_functions import *


# -----TEST CONVERTING VELOCITIES TO PM FUNCTION---------------------------------------------------------------------------------------------------------
# Galactic coordinate system is centered around sun
galactic_coord = SkyCoord(l=0*u.deg, b=90*u.deg, frame="galactic")
equatorial_coord = galactic_coord.icrs
eq_ra = equatorial_coord.ra.deg
eq_dec = equatorial_coord.dec.deg

# STAR MOVING PARALLEL TO SUN
def test_parallel_star():
    star4 = convert_vel_to_pm(ra=eq_ra, dec=eq_dec, distance=100, v_R=U0, v_phi=V0+Theta_0, v_z=W0)
    np.testing.assert_allclose(star4.pm_ra_cosdec.to_value(), 0, atol=1e-5)
    np.testing.assert_allclose(star4.pm_dec.to_value(), 0, atol=1e-5)
    np.testing.assert_allclose(star4.radial_velocity.to_value(), 0, atol=1e-7) 
    print("Star moving parallel to sun test passed")

# STAR AT EQUATOR WITH V_PHI = 0
def test_zero_v_phi_star():
    galactic_coord = SkyCoord(l=0*u.deg, b=0*u.deg, frame="galactic")
    equatorial_coord = galactic_coord.icrs
    eq_ra = equatorial_coord.ra.deg
    eq_dec = equatorial_coord.dec.deg
    star5 = convert_vel_to_pm(ra=eq_ra, dec=eq_dec, distance=100, v_R=U0, v_phi=0, v_z=W0)
    star5_gal = star5.transform_to(Galactic())
    np.testing.assert_allclose(star5_gal.pm_b.to_value(), 0, atol=1e-3)
    np.testing.assert_allclose(star5_gal.radial_velocity.to_value(), 0, atol=1e-3) 
    print("Star at equator with v_phi=0 test passed")


test_parallel_star()
test_zero_v_phi_star()