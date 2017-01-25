# Code to generate the geometry for my simulations
# Author: D. Michael Senter
#         Miller Lab, UNC-CH
# Date:   September 2016
from numpy import linspace
import re
from geo_obj import *

def findParams(fname):
# function that reads the grid parameters from the input2d file
# returns a dictionary
    params = {}
    with open(fname, 'r') as inFile:
        for line in inFile:
            if 'Nx' in line:
                for s in line.split():
                    if s.isdigit():
                        params['Nx'] = int(s)
            elif 'Lx' in line:
                for s in line.split():
                    if re.match("^\d+?(\.\d+)?$", s) is not None:
                        params['Lx'] = float(s)
    return params

def main():
    projname = "thrips"

    # the range we need
    params = findParams('input2d')
    Lx = params['Lx']
    Nx = params['Nx']
    xgrid = linspace(0,Lx,num=Nx)

    # create the lagragian points
    # starting with one mass, so only one point is needed
    lagpt1 = Vertex(0.5*(xgrid[15]+xgrid[16]),0.5*(xgrid[25]+xgrid[26]))
    vvec = [lagpt1]

    # make the mass point
    mpt1 = Mass(1,1e6,6e-8)
    mvec = [mpt1]

    # write the geometry:
    writeFile(projname, vvec)
    writeFile(projname, mvec)


if __name__=='__main__': main()
