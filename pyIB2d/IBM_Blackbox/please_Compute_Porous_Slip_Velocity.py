'''-------------------------------------------------------------------------

 IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based off of
	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

 Author: Nicholas A. Battista
 Email:  nick.battista@unc.edu
 Date Created: May 27th, 2015\
 Initial Python 3.5 port by: Christopher Strickland
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
from math import sqrt

################################################################################
#
# FUNCTION: computes Porous Slip Velocity based on Darcy's Law,
#
#           U_porous = -alpha <F_Lag,hat{n}> / | d vec{X}_Lag /ds |
#
################################################################################

def please_Compute_Porous_Slip_Velocity(ds,xLag,yLag,porous_info,F_Lag):
    ''' Computes porous slip velocity based on Darcy's Law
    
    Args:
        ds:
        xLag: vector of x-Pts associated w/ x-Lagrangian pts
        yLag: vector of y-Pts associated w/ y-Lagrangian pts
        porous_info:  col 1: lag-ids for porous media
                      col 2: x-Lag pts for lag-ids 
                      col 3: y-Lag pts for lag-ids
                      col 4: porosity coefficient
                      col 5: porous pt. stencil flag
        F_Lag:
        
    Returns:
        Por_Mat:
        nX:
        nY:'''

    # # of porous media pts.
    Np = porous_info.shape[0]

    # Initialize storage
    #Por_X = np.zeros(len(xLag))
    #Por_Y = np.zeros(len(xLag))

    # Compute Lagrangian Derivatives
    xL_s,yL_s = give_Me_Lagrangian_Derivatives(ds,Np,porous_info)

    # Compute Normal Vector (unit normals)
    nX,nY,sqrtNorm = give_Me_Lagrangian_Normal_Vectors(xL_s,yL_s)


    # Compute Porous Slip Velocity
    Up_X = - porous_info[:,3]*F_Lag[porous_info[:,0].astype('int'),0]*nX / sqrtNorm
    Up_Y = - porous_info[:,3]*F_Lag[porous_info[:,0].astype('int'),1]*nY / sqrtNorm
    Por_Mat = np.array([Up_X, Up_Y]).T

    # Store porous slip velocities in appropriate vector for adding to 
    #   current Lagrangian Velocity Computation
    #Por_X( porous_info[:,0] ) = Up_X
    #Por_Y( porous_info[:,0] ) = Up_Y
    #Por_Mat = np.array([Por_X, Por_Y])

    return (Por_Mat,nX,nY)

################################################################################
#
# FUNCTION: computes Lagrangian Derivatives
#
################################################################################

def give_Me_Lagrangian_Derivatives(ds,Np,porous_info):
    ''' Computes Lagrangian derivatives
    
    Args:
        ds:             # Lagrangian spacing
        Np:             # of Porous Pts.
        porous_info:    # Structure holding all porous information
               col 1:   # Porous IDs
               col 2:   # Porous x-Values
               col 3:   # Porous y-Values
               col 4:   # Porosity Coefficient
               col 5:   # Porous Pt. stencil flag
        
    Returns:
        xL_s:
        yL_s:'''

    xL = porous_info[:,1]    # x-Values for Porous points 
    yL = porous_info[:,2]    # y-Values for Porous points
    c = porous_info[:,4];    # choice of where it falls on stencil

    # Allocate storage for derivatives
    xL_s = np.zeros(Np)
    yL_s = np.zeros(Np)

    for i in range(0,Np):
        if c[i] == -2:
            xL_s[i] = ( -25/12*xL[i] + 4*xL[i+1] - 3*xL[i+2] + 4/3*xL[i+3] - 1/4*xL[i+4] ) / ds
            yL_s[i] = ( -25/12*yL[i] + 4*yL[i+1] - 3*yL[i+2] + 4/3*yL[i+3] - 1/4*yL[i+4] ) / ds
        
        elif c[i] == -1:
            xL_s[i] = ( -0.25*xL[i-1] - 5/6*xL[i] + 1.5*xL[i+1] - 0.5*xL[i+2] + 1/12*xL[i+3] ) / ds;
            yL_s[i] = ( -0.25*yL[i-1] - 5/6*yL[i] + 1.5*yL[i+1] - 0.5*yL[i+2] + 1/12*yL[i+3] ) / ds;

        elif c[i] == 0:
            xL_s[i] = ( 1/12*xL[i-2] - 2/3*xL[i-1] + 2/3*xL[i+1] - 1/12*xL[i+2] ) / ds;
            yL_s[i] = ( 1/12*yL[i-2] - 2/3*yL[i-1] + 2/3*yL[i+1] - 1/12*yL[i+2] ) / ds;
        
        elif c[i] == 1:
            xL_s[i] = ( -1/12*xL[i-3] + 0.5*xL[i-2] - 1.5*xL[i-1] + 5/6*xL[i] + 0.25*xL[i+1] ) / ds;
            yL_s[i] = ( -1/12*yL[i-3] + 0.5*yL[i-2] - 1.5*yL[i-1] + 5/6*yL[i] + 0.25*yL[i+1] ) / ds;
        
        elif c[i] == 2:
            xL_s[i] = ( 0.25*xL[i-4] - 4/3*xL[i-3] + 3*xL[i-2] - 4*xL[i-1] + 25/12*xL[i] ) / ds;
            yL_s[i] = ( 0.25*yL[i-4] - 4/3*yL[i-3] + 3*yL[i-2] - 4*yL[i-1] + 25/12*yL[i] ) / ds;


    # IF CLOSED POROUS STRUCTURE!!!
    #xL_s[0] = ( xL[1] - xL[-1] ) / (2*ds)
    #yL_s[0] = ( yL[1] - yL[-1] ) / (2*ds)
    #
    #xL_s[1:-1] = ( xL[2:] - xL[:-2] ) / (2*ds) 
    #yL_s[1:-1] = ( yL[2:] - yL[:-2] ) / (2*ds)
    #
    #xL_s[-1] = ( xL[0] - xL[-2] ) / (2*ds)
    #yL_s[-1] = ( yL[0] - yL[-2] ) / (2*ds)
    
    return (xL_s,yL_s)

################################################################################
# #
# # FUNCTION: compute Porous Connections
# #
# ##############################################################################

# def give_Me_Porous_Connections(ds,Np,porous_info):
    # '''Compute porous connections
    
    # Args:
        # ds:
        # Np:
        # porous_info:
        
    # Returns:
        # connects_Por:'''

    # xL = porous_info[:,1]
    # yL = porous_info[:,2]

    # start = 1
    # last = Np-1
    # prev = last
    # for ii in range(Np):
        
        # x0 = xL[ii]; y0 = yL[ii]      # central node value
        # xm1= xL[prev]; ym1 = yL[prev] # 'previous' values
        # xp1 =xL[ii+1]; yp1=yL[ii+1];  # 'next' values
        
        # distL = sqrt( (x0-xm1)**2 + (y0-ym1)**2 )
        # distR = sqrt( (x0-xp1)**2 + (y0-yp1)**2 )
        
        
    # return connects_Por

################################################################################
#
# FUNCTION: computes Lagrangian UNIT Normal Vectors
#
################################################################################

def give_Me_Lagrangian_Normal_Vectors(xL_s,yL_s):
    '''Computes Lagrangian UNIT normal vectors
    
    Args:
        xL_s:
        yL_s:
        
    Returns:
        nX:
        nY:
        sqrtN:'''

    sqrtN = np.sqrt( (xL_s)**2 + (yL_s)**2 )

    nX = ( yL_s ) / sqrtN
    nY = ( -xL_s) / sqrtN

    return (nX,nY,sqrtN)