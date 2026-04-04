"""Utilities for reading/writing Fortran binary files used by cfram_rrtmg."""

import numpy as np


def write_fortran_direct(filename, data, recl_bytes=None):
    """Write data to Fortran direct-access unformatted binary file.

    The cfram_rrtmg.f90 code uses direct-access I/O:
        open(unit=N, file=fname, access='direct', form='unformatted', recl=R)
        read(unit, rec=itime) data

    With gfortran -fdefault-real-8, real values are 8 bytes (float64).

    Args:
        filename: output file path
        data: numpy array (float64) to write as a single record
        recl_bytes: record length in bytes (if None, computed from data size)
    """
    arr = np.asarray(data, dtype=np.float64)
    arr.tofile(filename)


def read_fortran_direct(filename, shape, dtype=np.float64):
    """Read Fortran direct-access unformatted binary file.

    Args:
        filename: input file path
        shape: shape of the output array
        dtype: data type (default float64 for -fdefault-real-8)

    Returns:
        numpy array with given shape
    """
    data = np.fromfile(filename, dtype=dtype)
    return data.reshape(shape)


def write_scalar_direct(filename, value):
    """Write a single scalar to Fortran direct-access file."""
    np.array([value], dtype=np.float64).tofile(filename)


def write_3d_profile(filename, data_3d):
    """Write a 3D profile (nlev, nlat, nlon) in Fortran loop order.

    Fortran reads as: (((var(ilev,ilat,ilon), ilon=1,nlon), ilat=1,nlat), ilev=1,nlev)
    This is just column-major (Fortran) order, which for C-order arrays
    means we need to transpose and write.

    For the single-column case (nlat=1, nlon=1), this is just the 1D profile.
    """
    arr = np.asarray(data_3d, dtype=np.float64)
    # Fortran order: fastest varying index is rightmost in the read statement
    # read(unit, rec=itime) (((var(ilev,ilat,ilon), ilon=1,nlon), ilat=1,nlat), ilev=1,nlev)
    # This reads: for each ilev, for each ilat, for each ilon
    # In memory: contiguous by (ilon, ilat, ilev) = Fortran column-major of (nlev, nlat, nlon)
    arr.tofile(filename)


def write_4d_aerosol(filename, data_4d):
    """Write 4D aerosol data (nlev, nlat, nlon, nbnd) in Fortran loop order.

    Fortran reads as: ((((var(ilev,ilat,ilon,ib), ib=...), ilon=1,nlon), ilat=1,nlat), ilev=1,nlev)
    """
    arr = np.asarray(data_4d, dtype=np.float64)
    arr.tofile(filename)
