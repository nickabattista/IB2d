#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
#     fluid-structure interaction models. This version of the code is based off of
#    Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
# Email:  nick.battista@unc.edu
# Date Created: May 27th, 2015
# Institution: UNC-CH
#
# This code is capable of creating Lagrangian Structures using:
#     1. Springs
#     2. Beams (*torsional springs)
#     3. Target Points
#    4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
#
# One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
# 
# There are a number of built in Examples, mostly used for teaching purposes. 
# 
# If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
#
#--------------------------------------------------------------------------------------------------------------------#

import numpy as np
#import sys
#import os
from read_vtk_data import read_Eulerian_Data_From_vtk

#################################################################################
#
# FUNCTION: gives (x,y) positions of the immersed boundary at a single step
#          
################################################################################

def import_Eulerian_Data(path,numSim,Eulerian_Flags):
 
        #
        # EULERIAN FLAGS FOR WHAT GETS SPIT OUT #
        # 
        #Eulerian_Flags(0):   OMEGA
        #Eulerian_Flags(1):   PRESSURE
        #Eulerian_Flags(2):   uMAG
        #Eulerian_Flags(3):   uX (mag. x-component of velocity)
        #Eulerian_Flags(4):   uY (mag. x-component of velocity)
        #Eulerian_Flags(5):   uVEC (x,y-components of velocity: U,V)
        #Eulerian_Flags(6):   Fx (x-component of force )
        #Eulerian_Flags(7):   Fy (y-component of force)
        #

    # read in Vorticity #
    if Eulerian_Flags[0]:
        strChoice = 'Omega'; getxy = True
        Omega,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        Omega=[]

    # read in Pressure #
    if Eulerian_Flags[1]:
        strChoice = 'P'; getxy = True
        P,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        P = []

    # read in Velocity Magnitude #
    if Eulerian_Flags[2]:
        strChoice = 'uMag'; getxy = True
        uMag,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        uMag = []
 
    # read in x-directed Velocity Magnitude #
    if Eulerian_Flags[3]:
        strChoice = 'uX'; getxy = True
        uX,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        uX = []
 
    # read in y-directed Velocity Magnitude #
    if Eulerian_Flags[4]:
        strChoice = 'uY'; getxy = True
        uY,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        uY = []

    # read in x-directed Forces #
    if Eulerian_Flags[6]:
        strChoice = 'Fx'; getxy = True
        Fx,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        Fx = []

    # read in y-directed Forces #
    if Eulerian_Flags[7]:
        strChoice = 'Fy'; getxy = True
        Fy,x,y = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        Fy = []
 
    # read in Velocity Field #
    if Eulerian_Flags[5]:
        strChoice = 'u'; getxy = True
        U,V = read_Eulerian_Data_From_vtk(path,numSim,strChoice,getxy)
    else:
        U = []
        V = []

    # default for x,y values
    inds = np.array([0,1,2,3,4,6,7])
    if max(Eulerian_Flags[inds])==0:
        x=[]
        y=[]

    return x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy

