#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: rahul.mahajan@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z rahul.mahajan@noaa.gov $
###############################################################

__author__  = 'Barry Baker'
__email__   = 'barry.baker@noaa.gov'
__license__ = 'GPL'

'''
Simple utility to convert NCEP grib2 aerosol file into a netCDF4 file
Utilizes wgrib2 utility from NWPROD
'''

import os
from glob import glob
import sys
import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def execute_subprocess(cmd, verbose=False):
    '''
    call the subprocess command
    :param cmd: command to execute
    :param type: str
    '''
    if verbose:
        print( 'Executing: %s' % cmd)

    try:
        out = subprocess.check_output(cmd,shell=True)
    except subprocess.CalledProcessError as e:
        raise subprocess.CalledProcessError('',cmd,output=e.output)
    return


def get_exec_path(exec_name, verbose=False):
    '''
    get the full path to a given executable name
    :param exec_name: executable to fine
    :param type: str
    '''

    exec_path_def = '/nwprod2/grib_util.v1.0.0/exec/%s' % exec_name

    exec_path = find_executable(exec_name)
    if exec_path is None:
        exec_path = exec_path_def

    if verbose:
        print( '%s: %s' % (exec_name, exec_path))

    return exec_path

def chdir(fname):
#    dir_path = os.path.dirname(os.path.realpath(fname))
    dir_path = os.path.dirname(fname)
    os.chdir(dir_path)
    return os.path.basename(fname)


def change_file(finput,verbose=False):
    # first change directory and get file name
#    fname = chdir(finput)
    fname = finput
    # this will create 3 files and append to them
    wgrib2 = get_exec_path('wgrib2', verbose=verbose)

    # ENTIRE ATMOSPHERE GRIB LAYER
    #cmd = '%s %s -match "entire atmosphere:" -nc_nlev 1 -append -set_ext_name 1 -netcdf %s.entire_atm.nc' % (wgrib2, fname, fname)
    #execute_subprocess(cmd, verbose=verbose)
    # 1 hybrid level:
    #cmd = '%s %s -match "1 hybrid level:" -append -set_ext_name 1 -netcdf %s.hybrid.nc' % (wgrib2, fname, fname)
    #execute_subprocess(cmd, verbose=verbose)
    # surface:
    #cmd = '%s %s -match "surface:" -nc_nlev 1 -append -set_ext_name 1 -netcdf %s.surface.nc' % (wgrib2, fname, fname)
    #execute_subprocess(cmd, verbose=verbose)

    # VIIRS SURFACE GREENESS FRACTION
    cmd = '%s %s -match "surface:" -nc_nlev 1 -append -set_ext_name 1 -netcdf %s.nc' % (wgrib2, fname, fname)
    execute_subprocess(cmd, verbose=verbose)


if __name__ == '__main__':

    parser = ArgumentParser(description='convert grib2 file to netCDF4 file', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--files', help='input grib2 file name', type=str, required=True)
    parser.add_argument('-v', '--verbose', help='print debugging information', action='store_true', required=False)
    args = parser.parse_args()

    finput = args.files
    verbose = args.verbose

    files = sorted(glob(finput))
    for i,j in enumerate(files):
#        files[i] = os.path.realpath(j)
        files[i] = j
    if len(files) == 1:
        finput = files[0]
        change_file(finput,verbose=verbose)
    else:
        for i in files:
            finput = i
            print('FINPUT -> ',finput)
            change_file(finput,verbose=verbose)

    sys.exit(0)
