"""
utils
------

Useful tools used by the package.

.. autosummary::

   cart2sph
   sph2cart
   rot_xz
   rot_xz_sph

"""

import numpy as np

_DEG2RA = np.pi / 180

def cart2sph(vec):
    """ Converts cartesian [x,y,z] to spherical [r, theta, phi] coordinates 
    (in degrees).
    
    Parameters
    ----------
    vec: array
        x, y, z

    Returns
    -------
    array
        [r, theta, phi]
    """
    x, y ,z = vec
    v = np.sqrt(x**2 + y**2 + z**2)
    return np.array([v,
                    (np.arctan2(y,x) / _DEG2RA + 180) % 360 - 180, 
                     np.arcsin(z/v) / _DEG2RA])


def sph2cart(vec):
    """ Converts spherical coordinates [r, theta, phi]
    to cartesian coordinates [x,y,z].
    
    Parameters
    ----------
    vec: array
        r, theta, phi ; angles in degrees

    Returns
    -------
    array
        [x, y, z]
    """
    v, l, b = vec[0], np.asarray(vec[1])*_DEG2RA, np.asarray(vec[2])*_DEG2RA
    return np.asarray([v*np.cos(b)*np.cos(l), 
                       v*np.cos(b)*np.sin(l), 
                       v*np.sin(b)])  
     
def rot_xz(vec, theta):
    """ Rotates cartesian vector v [x,y,z] by angle theta around axis (0,1,0) 

    Parameters
    ----------
    vec: array
        x, y, z

    theta: float
        angle in degree

    Returns
    -------
    array
        rotated x, y, z
    """
    return [vec[0]*np.cos(theta*_DEG2RA) - vec[2]*np.sin(theta*_DEG2RA),
            vec[1][None,:],
            vec[2]*np.cos(theta*_DEG2RA) + vec[0]*np.sin(theta*_DEG2RA)]

def rot_xz_sph(l, b, theta):
    """ Rotate spherical coordinate (l,b = theta, phi) by angle theta around axis (0,1,0)
    (calls does to rot_xz and cart2sph)
    
    Parameters
    ----------
    l, b: float
       spherical coordinate

    theta: float
        angle in degree

    Returns
    -------
    array
        [r, theta, phi]
    """
    v_rot = rot_xz( sph2cart([1,l,b]), theta)
    return cart2sph(v_rot)[1:]
