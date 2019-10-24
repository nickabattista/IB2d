'''-------------------------------------------------------------------------

 IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based off of
	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

 Author: Nicholas A. Battista
 Email:  nick.battista@unc.edu
 Date Created: May 27th, 2015\
 Python 3.5 port by: Christopher Strickland
 Institution: UNC-CH

 This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")

 One is able to update those Lagrangian Structure parameters, e.g., 
 spring constants, resting lengths, etc
 
 There are a number of built in Examples, mostly used for teaching purposes. 
 
 If you would like us to add a specific muscle model, 
 please let Nick (nick.battista@unc.edu) know.

----------------------------------------------------------------------------'''

import numpy as np
import sys
from sys import platform
if platform == 'darwin': # OSX backend does not support blitting
    import matplotlib
    matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.rcParams['image.cmap'] = 'viridis'
from math import sqrt
import warnings

###########################################################################
#
# FUNCTION: Plots the Lagrangian structure with:
#           (1): the background velocity field
#           (2): the Lagrangian structure itself
#           (3): the vorticity field in a colormap
#
###########################################################################

def please_Plot_Results(ds,X,Y,U,V,vort,uMag,p,chiX,chiY,\
    lagPlot,velPlot,vortPlot,pressPlot,uMagPlot,firstPrint,loc,diffy,spacing):
    '''Plots the Lagrangian structure w/ velocity field and vorticity
    
    Args:
        ds:
        X,Y:  (x,y) values
        U,V:  x-directed, y-directed velocities respectively
        vort: vorticity
        uMag: magnitude of velocity
        p:    pressure
        chiX,chiY:
        lagPlot:   - if you want lag. point ONLY plot = 1
        velPlot:   - if you want velocity plot = 1
        vortPlot:  - if you want vorticity plot = 1
        pressPlot: - if you want pressure plot = 1
        uMagPlot:  - if you want mag. velocity plot = 1
        firstPrint:  flag, 1 or 0
        loc:  Altered internally!
        diffy:   Altered internally!
        spacing:
        
    Returns:
        loc:
        diffy:'''

    #
    # Assumption: Assuming chiX and chiY are column vectors
    # Assumption: Assuming chiX(i+1)-chiX(i) < .5 and chiY(i+1)-chiY(i) < .5, for all points that don't cross the boundary
    #
    # This script has been adapted from the Mat-IB code by Hiens and Stockie, 2014.
    #

    Lx = X[0,-1]+X[0,1]
    Ly = Y[-1,0]+Y[1,0]

    #Shift so inside of interval [0,Lx] or [0,Ly]
    chiX = chiX % Lx
    chiY = chiY % Ly

    #Connect Geometry Appropriately At Beginning Of Simulation then Store it!
    if firstPrint:
        diffX = np.abs(np.roll(chiX,-1)-chiX)
        diffY = np.abs(np.roll(chiY,-1)-chiY)
        
        # assuming x value can't change more than Lx/2 (without crossing boundary)
        locX = np.argwhere(diffX > spacing).flatten()
        # assuming y value can't change more than Ly/2 (without crossing boundary)
        locY = np.argwhere(diffY > spacing).flatten()
        loc = np.unique(np.hstack((locX,locY)))+1 #automatically sorts the output
        #adds one so that index points to after the jump. This works better for
        #   Python 0-base indexing, with endpt not included

        diffy = sqrt( (chiX[-1]-chiX[0] )**2 + ( chiY[-1]-chiY[0] )**2 )

    plt.clf() #Clear previous plots :)




    plt.figure(1) 
    numPlots = lagPlot+velPlot+vortPlot+pressPlot+uMagPlot

    ct = 1

    # # # # # PLOTS LAGRANGIAN POINTS ONLY (if selected) # # # # #

    if lagPlot:
        plt.subplot(1,numPlots,ct)
        plt.axis([0, Lx, 0, Ly])
        plt.title('LAGRANGIAN PTS')
        plt.xlabel('x'); plt.ylabel('y')

        loc = np.hstack((0,loc,chiX.size-1))
        for ii in range(1,loc.size):     #=2:length(loc)
            plt.plot(chiX[loc[ii-1]:loc[ii]],chiY[loc[ii-1]:loc[ii]],'m',linewidth=3)
        if diffy < 5*ds:
            xTemp = [chiX[0], chiX[-1]]
            yTemp = [chiY[0], chiY[-1]]
            plt.plot(xTemp[0:2],yTemp[0:2],'m',linewidth=3)

        plt.axis('square')
        plt.axis([0, Lx, 0, Ly])

        ct+=1

    # # # # # PLOTS VORTICITY + LAG Pts. (if selected) # # # # #

    if vortPlot:
        
        plt.subplot(1,numPlots,ct)
        #
        plt.axis([0, Lx, 0, Ly])
        plt.title('VORTICITY')
        plt.xlabel('x'); plt.ylabel('y')

        #Compute Vorticity and Plot It against Lagrangian Grid!
        x = X[0,:]; y = Y[:,0]
        plt.contourf(x,y,np.flipud(np.rot90(vort)),10)

        loc = np.hstack((0,loc,chiX.size))
        for ii in range(1,loc.size):
            plt.plot(chiX[loc[ii-1]:loc[ii]],chiY[loc[ii-1]:loc[ii]],'m',linewidth=3)
        if diffy < 5*ds:
            xTemp = [chiX[0], chiX[-1]]
            yTemp = [chiY[0], chiY[-1]]
            plt.plot(xTemp[0:2],yTemp[0:2],'m',linewidth=3)
        
        plt.axis('square')
        plt.axis([0, Lx, 0, Ly])


        ct+=1

    # # # # # PLOTS PRESSURE + LAG Pts. (if selected) # # # # #

    if pressPlot:
        
        plt.subplot(1,numPlots,ct)
        #
        plt.axis([0, Lx, 0, Ly])
        plt.title('PRESSURE')
        plt.xlabel('x'); plt.ylabel('y') 

        #Use Pressure and Plot It against Lagrangian Grid!
        x = X[0,:]; y = Y[:,0]
        plt.contourf(x,y,p,6)

        loc = np.hstack((0,loc,chiX.size))
        for ii in range(1,loc.size):
            plt.plot(chiX[loc[ii-1]:loc[ii]],chiY[loc[ii-1]:loc[ii]],'m',linewidth=3)
        if diffy < 5*ds:
            xTemp = [chiX[0], chiX[-1]]
            yTemp = [chiY[0], chiY[-1]]
            plt.plot(xTemp[0:2],yTemp[0:2],'m',linewidth=3)

        plt.axis('square')
        plt.axis([0, Lx, 0, Ly])


        ct+=1

    # # # # # PLOTS MAGNITUDE OF VELOCITY + LAG Pts. (if selected) # # # # #

    if uMagPlot:
        
        plt.subplot(1,numPlots,ct)
        #
        plt.axis([0, Lx, 0, Ly])
        plt.title('MAGNITUDE OF VELOCITY')
        plt.xlabel('x'); plt.ylabel('y');

        #Use Mag. Velocity and Plot It against Lagrangian Grid!
        x = X[0,:]; y = Y[:,0]
        plt.contourf(x,y,uMag,6)

        loc = np.hstack((0,loc,chiX.size))
        for ii in range(1,loc.size):
            plt.plot(chiX[loc[ii-1]:loc[ii]],chiY[loc[ii-1]:loc[ii]],'m',linewidth=3)
        if diffy < 5*ds:
            xTemp = [chiX[0], chiX[-1]]
            yTemp = [chiY[0], chiY[-1]]
            plt.plot(xTemp[0:2],yTemp[0:2],'m',linewidth=3)

        plt.axis('square')
        plt.axis([0, Lx, 0, Ly])


        ct+=1

    # # # # # PLOTS VELOCITY FIELD + LAG Pts. (if selected) # # # # #

    if velPlot:
        
        plt.subplot(1,numPlots,ct)
        #
        plt.axis([0, Lx, 0, Ly])
        plt.title('VELOCITY')
        plt.xlabel('x'); plt.ylabel('y')

        plt.quiver(X,Y,U,V) #Print Velocity Field

        loc = np.hstack((0,loc,chiX.size))
        for ii in range(1,loc.size):
            plt.plot(chiX[loc[ii-1]:loc[ii]],chiY[loc[ii-1]:loc[ii]],'m',linewidth=3)
        if diffy < 5*ds:
            xTemp = [chiX[0], chiX[-1]]
            yTemp = [chiY[0], chiY[-1]]
            plt.plot(xTemp[0:2],yTemp[0:2],'m',linewidth=3)

        plt.axis('square')
        plt.axis([0, Lx, 0, Ly])

        
        #ct+=1
    
    plt.box(on=True)
    
    plt.draw()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.pause(0.0001) #no idea why this is necessary, but it is

    

    return (loc, diffy)
