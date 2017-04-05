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
#           (1) This code analyzes viz_IB2d code for channel flow in
#               /data_analysis/Example_For_Data_Analysis/Example_Flow_In_Channel/viz_IB2d
#           (2) Produces a plot of cross-sectional mag. of velocity for
#               different points along the channel at three times.
#           (3) USER-DEFINED functions are functions that users should make to
#               analyze their specific data sets
#           (4) MUST make sure to 'addpath' to where DA_Blackbox is, i.e.,
#               line 61
#           (5) MUST make sure to set path to desired dataset, i.e., in line
#               56
#          
#
################################################################################

def Example_Channel_Flow_Analysis():
    
    # TEMPORAL INFO FROM input2d #
    dt = 1e-4      # Time-step
    Tfinal = 0.015 # Final time in simulation
    pDump=50       # Note: 'print_dump' should match from input2d

    # DATA ANALYSIS INFO #
    start=1                                            # 1ST interval # included in data analysis
    finish=3                                           # LAST interval # included in data analysis
    dump_Times = np.array(range(start,finish+1))*pDump # Time vector when data was printed in analysis
        
    # SET PATH TO DESIRED viz_IB2d DATA %
    pathViz = '/Users/nick_battista/Desktop/IB2d/data_analysis/analysis_in_matlab/Example_For_Data_Analysis/Example_Flow_In_Channel/viz_IB2d/'
    pathForce = '/Users/nick_battista/Desktop/IB2d/data_analysis/analysis_in_matlab/Example_For_Data_Analysis/Example_Flow_In_Channel/hier_IB2d_data'

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


        #                                                                    #
        # *** USER DEFINED FUNCTIONS TO GET DESIRED ANALYSIS PT. INDICES *** #
        #                                                                    #
        if i == start:
            xPts = np.array([0.125,0.225,0.325,0.425])
            yPts = np.array([0.405,0.595])
            xInds,yInds = give_Desired_Analysis_Points(x,y,xPts,yPts)
            vel_data = np.zeros( (finish,len(yInds),len(xInds)) )

        #                                                                    #
        # ***** USER DEFINED FUNCTION TO SAVE DESIRED VELOCITY DATA *****    #
        #                                                                    #
        vel_data = store_Desired_Magnitude_Velocity_Data(uMag,vel_data,xInds,yInds,i)

    yVals = y[yInds]
    plot_Desired_Data(yVals,vel_data)

#################################################################################
#
# USER-FUNCTION: Finds desired analysis point indices
#
#        INPUTS: x: x-grid pts (row vector)
#                y: y-grid pts (column vector)
#                xPts: desired x-pts to analyze
#                yPts: desired y-pts to analyze
#
#################################################################################

def give_Desired_Analysis_Points(x,y,xPts,yPts):

    # Get x-Indices
    xInds = np.zeros(len(xPts))
    for i in range(0,len(xPts)):
        xPt = xPts[i]   # Get x-Pts
        k = 0           # Initialize loop iteration variable
        while x[k] < xPt:
            xInds[i] = k
            k = k+1

    # Get y-Indices
    yIndsAux = np.zeros(len(yPts))
    for i in range(0,len(yPts)):
        yPt = yPts[i]   # Get x-Pts
        k = 0           # Initialize loop iteration variable
        while y[k] < yPt:
            yIndsAux[i] = k
            k = k+1

    yInds = np.array(range(int(yIndsAux[0]),int(yIndsAux[1])))

    return xInds, yInds

#################################################################################
#
# USER-FUNCTION: Stores desired magnitude of velocity data
#
#        INPUTS: uMag: magnitude of velocity from simulation
#                vel_data: stored mag. of velocity data in 3D-matrix
#                xInds/yInds: indices of where to save data
#                i: ith time-step to store data from
#
#################################################################################

def store_Desired_Magnitude_Velocity_Data(uMag,vel_data,xInds,yInds,i):

# NOTE: vel_data(k,i): : k: magnitude of velocity data  (row)
#                      : j: xPt you're storing data for (column)
#                      : i: ith level you're storing in (level)

    for j in range(0,len(xInds)):
        for k in range(0,len(yInds)):
            vel_data[i-1,k,j] = uMag[yInds[k],xInds[j]]

    return vel_data



#################################################################################
#
# USER-FUNCTION: Plots magnitude of velocity data
#
#        INPUTS: vel_data: stored mag. of velocity data in 3D-matrix
#
#################################################################################

def plot_Desired_Data(yVals,vel_data):

    # Set Figure Size
    FigHand = plt.figure(1)
    #set(FigHand,'Position',[100,100,1024,895])
    #
    # Make Figure!
    #
    plt.subplot(3,1,1)
    mat = vel_data[0,:,:]
    maxVal = mat.max()
    plt.plot(yVals,vel_data[0,:,0],'*-');  plt.hold(True)
    plt.plot(yVals,vel_data[0,:,1],'r*-'); plt.hold(True)
    plt.plot(yVals,vel_data[0,:,2],'g*-'); plt.hold(True)
    plt.plot(yVals,vel_data[0,:,3],'k*-'); plt.hold(True)
    plt.axis([0.4,0.6,0,1.1*maxVal])
    #leg=plt.legend('x=0.125','x=0.175','x=0.225','x=0.275')
    plt.title('t=0.005')
    plt.ylabel('Mag. Velocity')
    plt.xlabel('y')
    #
    #
    plt.subplot(3,1,2)
    mat = vel_data[1,:,:]
    maxVal = mat.max()
    plt.plot(yVals,vel_data[1,:,0],'*-');  plt.hold(True)
    plt.plot(yVals,vel_data[1,:,1],'r*-'); plt.hold(True)
    plt.plot(yVals,vel_data[1,:,2],'g*-'); plt.hold(True)
    plt.plot(yVals,vel_data[1,:,3],'k*-'); plt.hold(True)
    plt.axis([0.4,0.6,0,1.1*maxVal])
    #leg=plt.legend('x=0.125','x=0.175','x=0.225','x=0.275')
    plt.title('t=0.01')
    plt.ylabel('Mag. Velocity')
    plt.xlabel('y')
    #
    #
    plt.subplot(3,1,3)
    mat = vel_data[2,:,:]
    maxVal = mat.max()
    plt.plot(yVals,vel_data[2,:,0],'*-');  plt.hold(True)
    plt.plot(yVals,vel_data[2,:,1],'r*-'); plt.hold(True)
    plt.plot(yVals,vel_data[2,:,2],'g*-'); plt.hold(True)
    plt.plot(yVals,vel_data[2,:,3],'k*-'); plt.hold(True)
    plt.axis([0.4,0.6,0,1.1*maxVal])
    #leg=plt.legend('x=0.125','x=0.175','x=0.225','x=0.275')
    plt.title('t=0.015')
    plt.ylabel('Mag. Velocity')
    plt.xlabel('y')


    plt.hold(False)
    plt.box(on=True)
    plt.draw()
    plt.pause(0.0001)


#################################################################################


if __name__ == "__main__":
    Example_Channel_Flow_Analysis()