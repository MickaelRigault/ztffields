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
import h5py

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

def ccdid_qid_to_rcid(ccdid, qid):
    """ provides the rcid corresponding to the ccdid and qid """
    return np.asarray(4*(ccdid - 1) + qid - 1, dtype="int")

def rcid_to_ccdid_qid(rcid):
    """ gets the ccdid and qid for the given rcid """
    qid = (rcid%4)+1
    ccdid  = ((rcid-(qid - 1))/4 +1)
    return np.asarray([ccdid, qid], dtype="int")

class ZTFCoords():
    _CCDSHAPE = (3080, 3072)
    _QUADRANTSHAPE = (..., ...)

    def __init__(self, coefs=None, h5file=None):
        """
        Initialize ZTFCoords.
        
        Either:
            - coefs: dict with keys 'u','v','x','y' (arrays shape [n_rcid,6])
        Or:
            - h5file: path to HDF5 file storing the same keys
        """
        if coefs is not None:
            self.update_polynomials(coefs)
        elif h5file is not None:
            self._load_polynomials_from_h5(h5file)
        else:
            raise ValueError("Must provide either `coefs` or `h5file`")

    def _load_polynomials_from_h5(self, h5file):
        """Load polynomial coefficients from HDF5 file"""
        with h5py.File(h5file, "r") as f:
            coefs = {key: np.array(f[key]) for key in ("u","v","x","y")}
        self.update_polynomials(coefs)

    def update_polynomials(self, coefs):
        """
        Validate and construct polynomial coefficients.
        Keys: 'u', 'v', 'x', 'y'
        Each value: array of shape (n_rcid, 6)
        """
        try:
            cu = np.asarray(coefs["u"])
            cv = np.asarray(coefs["v"])
            cx = np.asarray(coefs["x"])
            cy = np.asarray(coefs["y"])
        except KeyError as e:
            raise KeyError(f"Missing polynomial key: {e.args[0]}")

        for name, arr in zip(("u", "v", "x", "y"), (cu, cv, cx, cy)):
            if arr.ndim != 2 or arr.shape[1] != 6:
                raise ValueError(f"coefs['{name}'] must have shape (n_rcid, 6)")

        n_rcid = cu.shape[0]
        if not (cv.shape[0] == cx.shape[0] == cy.shape[0] == n_rcid):
            raise ValueError("Inconsistent rcid dimension across polynomials")

        self._coefs_u = cu
        self._coefs_v = cv
        self._coefs_x = cx
        self._coefs_y = cy

    def polynom2D(self, x, y, coeffs):
        """
        2D polynomial with 6 parameters
        """
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        coeffs = np.atleast_2d(coeffs)

        return (
            coeffs[:, 0]
            + coeffs[:, 1]*x
            + coeffs[:, 2]*y
            + coeffs[:, 3]*x*y
            + coeffs[:, 4]*x**2
            + coeffs[:, 5]*y**2
        )

    # ij -> xy -> uv
    def ij_to_xy(self, i, j, qid):
        """ Converts ccd coordinates (i,j) into quadrant coordinates (x,y) """
        # recall:
        # q2 | q1
        # --------
        # q3 | q4

        x = i
        y = j

        if qid in [1, 4]:
            x = i-self._CCDSHAPE[1]

        if qid in [1, 2]:
            y = j-self._CCDSHAPE[0]

        return np.stack([x, y])
        
    def xy_to_uv(self, x, y, rcid):
        """ Converts quadrant coordinates (x,y) into focal plane coordinates (u,v) """
        x = np.asarray(x)
        y = np.asarray(y)
        rcid = np.asarray(rcid)
        
        coeffs_u = self._coefs_u[rcid, :]
        coeffs_v = self._coefs_v[rcid, :]

        u = self.polynom2D(x, y, coeffs_u)
        v = self.polynom2D(x, y, coeffs_v)
        
        return u, v
    
    # uv -> xy -> ij
    def uv_to_xy(self, u, v, rcid):
        """ Converts focal plane coordinates (u,v) into quadrant coordinates (x,y) """
        u = np.asarray(u)
        v = np.asarray(v)
        rcid = np.asarray(rcid)

        coeffs_x = self._coefs_x[rcid, :]
        coeffs_y = self._coefs_y[rcid, :]

        x = self.polynom2D(u, v, coeffs_x)
        y = self.polynom2D(u, v, coeffs_y)
        
        return x, y
    
    def xy_to_ij(self, x, y, qid):
        """ Converts quadrant coordinates (x,y) into ccd coordinates (i,j) """
        # recall:
        # q2 | q1
        # --------
        # q3 | q4

        i = x
        j = y

        if qid in [1, 4]:
            i = x+self._CCDSHAPE[1]

        if qid in [1, 2]:
            j = y+self._CCDSHAPE[0]
        
        return np.stack([i, j])

    # Byproduct functions
    def uv_to_ij(self, u, v, ccdid):
        """ Converts (u,v) to (i,j) by selecting the correct quadrant """

        u = np.asarray(u)
        v = np.asarray(v)

        for qid in [1, 2, 3, 4]:
            rcid = ccdid_qid_to_rcid(ccdid, qid)

            x, y = self.uv_to_xy(u, v, rcid)

            # validity mask
            valid = (
                (x >= 0) & (x <= self._CCDSHAPE[1]) &
                (y >= 0) & (y <= self._CCDSHAPE[0])
            )

            if np.all(valid):
                return self.xy_to_ij(x, y, qid)

        raise ValueError("No valid quadrant found for given (u,v)")
    
    def ij_to_qid(self, i, j):
        if i >= self._CCDSHAPE[1]:   # right half
            if j >= self._CCDSHAPE[0]:
                return 1  # top-right
            else:
                return 4  # bottom-right
        else:                         # left half
            if j >= self._CCDSHAPE[0]:
                return 2  # top-left
            else:
                return 3  # bottom-left

    def ij_to_uv(self, i, j, ccdid):
        """ Converts ccd coordinates (i,j) into focal plane coordinates (u, v) """
        qid = self.ij_to_qid(i, j)
        x, y = self.ij_to_xy(i,j, qid=qid)
        return self.xy_to_uv(x, y, ccdid_qid_to_rcid(ccdid,qid))
