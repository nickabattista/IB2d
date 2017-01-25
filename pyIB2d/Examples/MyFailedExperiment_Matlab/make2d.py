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

def WingTips():
    # function to find wing tip coords based on wing length, angle with horiz.
    # and position of center point
    pass

def WingY():
    # return y-val of an x on the wing
    pass

def CreateWing():
    # given wing tips and center mess as well as Ds, generate wing vertices
    pass

def main():
    projname = "thrips"

    # the range we need
    params = findParams('input2d')
    Lx = params['Lx']
    Nx = params['Nx']
    Ds = 0.5*Lx
    y = 0.7
    xgrid = linspace(0.1,0.6, num=12)
    veclist = []

    for x in xgrid:
        point = Vertex(x,y)
        veclist.append(point)

    splist = []
    for i in range(len(veclist)-1):
        spr = Spring(i+1,i+2)
        splist.append(spr)

    massvec = []
    for i in range(len(veclist)):
        masspt = Mass(i,1e6,1e-3)

    # write the geometry:
    writeFile(projname, veclist)
    writeFile(projname, splist)
    #writeFile(projname, mvec)


if __name__=='__main__': main()
