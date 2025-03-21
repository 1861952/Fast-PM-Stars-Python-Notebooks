import numpy as np
import pandas as pd
from scipy.stats import norm, powerlaw
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import Galactocentric, galactocentric_frame_defaults, ICRS, SkyCoord, Galactic

# Necessary Constants
# km/s
U0 = 11.1 
V0 = 12.24
W0 = 7.25
Theta_0 = 220
galcen_params = dict(galcen_v_sun=np.array([U0, V0+Theta_0, W0]) * u.km/u.s, galcen_distance=8*u.kpc)

# ra and dec given in equatorial coords. Velocities in galactocentric coords in cylindrical system
def convert_vel_to_pm(ra, dec, distance, v_R, v_phi, v_z):
    star_eq_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, distance=distance*u.pc, frame=ICRS())
    star_galcen_coord = star_eq_coord.transform_to(Galactocentric(**galcen_params))
    phi = star_galcen_coord.cylindrical.phi
    
    v_x = v_R*np.cos(np.pi*u.rad - phi) - v_phi*np.sin(np.pi*u.rad - phi)
    v_y = v_R*np.sin(np.pi*u.rad - phi) + v_phi*np.cos(np.pi*u.rad - phi)
    #  v_z stays the same
    
    star_galcen_coord_vel = Galactocentric(x=star_galcen_coord.x, y=star_galcen_coord.y, z=star_galcen_coord.z, v_x=v_x*u.km/u.s, v_y=v_y*u.km/u.s, v_z=v_z*u.km/u.s, **galcen_params)
    star_eq_coord_vel = star_galcen_coord_vel.transform_to(ICRS())
    return star_eq_coord_vel


# Cartesian to Cylindrical Velocities
def v_xyz_to_cyl(x, y, z, v_x, v_y, v_z):
    v_r = (x*v_x + y*v_y) / np.sqrt(x**2 + y**2)
    v_phi = (-y*v_x + x*v_y) / np.sqrt(x**2 + y**2)
    print("Sun's cylindrical velocities: ", v_r, v_phi, v_z)
    return v_r, v_phi, v_z

# Cylindrical to Cartesian Velocities
def v_cyl_to_xyz(v_r, v_phi, v_z, phi):
    v_x = -v_r*np.cos(phi) - v_phi*np.sin(phi)
    v_y = v_r*np.sin(phi) - v_phi*np.cos(phi)
    print("Sun's xyz velocities: ", v_x, v_y, v_z)
    return v_x, v_y, v_z
