#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled
# 	fluid-structure interaction models. This version of the code is based off of
#	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
#
# Author: Nicholas A. Battista
# Email:  nick.battista@unc.edu
# Date Created: May 27th, 2015
# Institution: UNC-CH
#
# This code is capable of creating Lagrangian Structures using:
# 	1. Springs
# 	2. Beams (*torsional springs)
# 	3. Target Points
#	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
#
# One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
#
# There are a number of built in Examples, mostly used for teaching purposes.
#
# If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
#
#--------------------------------------------------------------------------------------------------------------------#

import numpy as np
import sys
import matplotlib.pyplot as plt

#################################################################################
#
# FUNCTION: This is where you can perform data analysis of a simulation
#           from your stored viz_IB2d information
#
#     Note:
#
#           (1) USER-DEFINED functions are functions that users should make to
#               analyze their specific data sets
#           (2) MUST make sure to 'addpath' to where DA_Blackbox is, i.e.,
#               line 67
#           (3) MUST make sure to set path to desired dataset, i.e., in line
#               63,64
#           (4) USUALLY EASIEST TO SIMPLY PUT viz_IB2d and hier_IB2d_data folders
#               in this directory
#           (5) MUST have VTK library installed, if using Python3 with Anaconda, 
#               please type <conda install -c menpo vtk> into your terminal
#
################################################################################

def Standard_Data_Analysis_Script():
    
    # TEMPORAL INFO FROM input2d #
    dt = 1e-4      # Time-step
    Tfinal = 0.015 # Final time in simulation
    pDump=50       # Note: 'print_dump' should match from input2d
    
    # DATA ANALYSIS INFO #
    start=1                                            # 1ST interval # included in data analysis
    finish=3                                           # LAST interval # included in data analysis
    dump_Times = np.array(range(start,finish+1))*pDump # Time vector when data was printed in analysis
    
    # SET PATH TO DESIRED viz_IB2d DATA %
    pathViz = 'viz_IB2d/'
    pathForce = 'hier_IB2d_data'
    
    # SET PATH TO DA_BLACKBOX %
    sys.path.append('../DA_Blackbox/')
    
    from give_Lag_Positions import give_Lag_Positions
    from import_Eulerian_Data import import_Eulerian_Data
    from import_Lagrangian_Force_Data import import_Lagrangian_Force_Data
    
    for i in range(start,finish+1):
        
        # Points to desired data viz_IB2d data file
        if i<10:
            numSim = '000'+str(i)
        elif i<100:
            numSim = '00'+str(i)
        elif i<1000:
            numSim = '0'+str(i)
        else:
            numSim = str(i)
    
    
        # Imports immersed boundary positions %
        xLag,yLag = give_Lag_Positions(pathViz,numSim)
        
        # Imports (x,y) grid values and ALL Eulerian Data %
        #                      DEFINITIONS
        #          x: x-grid                y: y-grid
        #       Omega: vorticity           P: pressure
        #    uMag: mag. of velocity
        #    uX: mag. of x-Velocity   uY: mag. of y-Velocity
        #    U: x-directed velocity   V: y-directed velocity
        #    Fx: x-directed Force     Fy: y-directed Force
        #
        #  Note: U(j,i): j-corresponds to y-index, i to the x-index
        #
        x,y,Omega,P,uMag,uX,uY,U,V,Fx,Fy = import_Eulerian_Data(pathViz,numSim)


        # Imports Lagrangian Pt. FORCE (magnitude) DATA %
        #                      DEFINITIONS
        #
        #      fX_Lag: forces in x-direction on boundary
        #      fY_Lag: forces in y-direction on boundary
        #       fLagMag: magnitude of force at boundary
        #   fLagNorm: magnitude of NORMAL force at boundary
        #   fLagTan: magnitude of TANGENT force at boundary
        #
    fX_Lag,fY_Lag,fLagMag,fLagNorm,fLagTan = import_Lagrangian_Force_Data(pathForce,numSim)

    
if __name__ == "__main__":
    Standard_Data_Analysis_Script()