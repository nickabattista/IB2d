#-------------------------------------------------------------------------------------------------------------------#
#
# IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
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
#   4. Mass Points
#   5. Porous Points
#	6. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
#   7. 3-Element Hill Muscle Model
#
# One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting-lengths, etc
# 
# There are a number of built in Examples, mostly used for teaching purposes. 
# 
# If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
#
#--------------------------------------------------------------------------------------------------------------------#

import numpy as np
import math

##################################################################################
#
# FUNCTION: computes normal and tangential forces on Lagrangian Boundary
#
#################################################################################

def please_Compute_Normal_Tangential_Forces_On_Lag_Pts(lagPts,F_Lag):
   
    # Computes Normal and Tangential Forces on Lagrangian Fiber Structure
    #Args:
    #    lagPts:    Nx2 matrix containing Lag. Pts (x-col0, y-col1)
    #    F_Lag:     Nx2 matrix containing Lag. Pt Forces (Fx-col0, Fy-col1)     
    #Returns:
    #    F_Tan_Mag F_Normal_Mag

    X = lagPts[:,0]             # Stores xLagrangian Pts 
    Y = lagPts[:,1]             # Stores yLagrangian Pts
    Npts = len(X);              # Computes # of Lagrangian Pts

    # Compute Normal/Tangential Vectors
    nX,nY,sqrtN = give_Me_Lagrangian_Normal_Vectors(Npts,X,Y)
    [tX,tY,sqrtN] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN);

    # Project Force Data onto Normal / Tangent Vectors
    [F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,F_Lag[:,0],F_Lag[:,1],nX,nY,tX,tY);
    
    # Compute Colormap Force Magnitude Scalings
    [F_Tan_Mag,F_Normal_Mag] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal);

    return (F_Tan_Mag, F_Normal_Mag)

#####################################################################################
#
# FUNCTION: computes Lagrangian UNIT Normal Vectors
#
#####################################################################################

def give_Me_Lagrangian_Normal_Vectors(Npts,X,Y):

    # Gives normal vectors from Lag. Structure at each Node
    # Args:
    #     Npts:     # of Lagrangian Pts (scalar)
    #     X,Y:      Nx1 matriices containing Lag. Pts, X and Y, respectively     
    # Returns:
    #     nX, nY, sqrtN

    # Compute Lagrangian Spacing
    ds = np.sqrt( ( X[3]-X[4] )*( X[3]-X[4] ) + ( Y[3]-Y[4] )*( Y[3]-Y[4] ) )

    # Gives Lagrangian Derivatives
    xL_s, yL_s = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y);

    sqrtN = np.sqrt( xL_s*xL_s + yL_s*yL_s )

    nX = ( yL_s ) / sqrtN;
    nY = ( -xL_s) / sqrtN;

    return (nX, nY, sqrtN)

#####################################################################################
#
# FUNCTION: computes Lagrangian Derivatives for Normal/Tangential Vector Computation
#
#####################################################################################

def give_Me_Lagrangian_Derivatives(ds,Npts,X,Y):

    # computes Lagrangian Derivatives for Normal/Tangential Computation
    # Args:
    #     ds:       Lagrangian grid spacing
    #     Npts:     # of Lagrangian Pts (scalar)
    #     X,Y:      Nx1 matriices containing Lag. Pts, X and Y, respectively     
    # Returns:
    #     xL_s, yL_s <- the X,Y derivative values

    xL_s = np.zeros( Npts ) # allocate storage for derivative values
    yL_s = np.zeros( Npts ) # allocate storage for derivative values

    for i in range(0,Npts):     # to iterate between 0,...,Npts-1
        if (i==0):
            xL_s[i] = ( X[i+1] - X[-1] ) / (2*ds)
            yL_s[i] = ( Y[i+1] - Y[-1] ) / (2*ds)
        elif i < Npts-1:
            xL_s[i] = ( X[i+1] - X[i-1] ) / (2*ds) 
            yL_s[i] = ( Y[i+1] - Y[i-1] ) / (2*ds)
        else:
            xL_s[i] = ( X[0] - X[-1-1] ) / (2*ds) 
            yL_s[i] = ( Y[0] - Y[-1-1] ) / (2*ds)

    return (xL_s, yL_s)



#####################################################################################
#
# FUNCTION: computes Lagrangian UNIT Tangent Vectors
#
#####################################################################################

def give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN):

    # Allocate storage
    tX = np.zeros( (Npts,1) );
    tY = np.zeros( (Npts,1) );

    # Rotate normal vectors to get tangent vectors
    ang = -math.pi/2; # Rotate CW by 90 degrees
    for i in range(0,Npts):
        tX[i] = nX[i]*math.cos(ang) - nY[i]*math.sin(ang);
        tY[i] = nX[i]*math.sin(ang) + nY[i]*math.cos(ang);

    return (tX , tY, sqrtN)
    
    # Testing
    #for i in range(0,Npts):
    #   test[i] = ( nX[i)*tX[i] + nY[i]*tY[i] ); 
    #   test2[i] = sqrt( tX[i]*tX[i] + tY[i]*tY[i] );
    #end


#####################################################################################
#
# FUNCTION: computes force vector projections onto the tangent and normal
#           vectors! 
#
#####################################################################################
    
def give_Tangent_and_Normal_Force_Projections(Npts,Fx,Fy,nX,nY,tX,tY):

    # Allocate Storage
    F_Tan = np.zeros( (Npts,2) )
    F_Normal = np.zeros( (Npts,2) )

    for i in range(0,Npts):

        # Compute dot product between force vector and tangent vector
        tanVec_dotProd = ( Fx[i]*tX[i] + Fy[i]*tY[i] ) / np.sqrt( tX[i]*tX[i] + tY[i]*tY[i] )
        F_Tan[i,0] = tanVec_dotProd * ( tX[i] )
        F_Tan[i,1] = tanVec_dotProd * ( tY[i] )
    
        # Compute dot product between force vector and normal vector
        normalVec_dotProd = ( Fx[i]*nX[i] + Fy[i]*nY[i] ) / np.sqrt( nX[i]*nX[i] + nY[i]*nY[i] )
        F_Normal[i,0] = normalVec_dotProd * ( nX[i] )
        F_Normal[i,1] = normalVec_dotProd * ( nY[i] )
    
    return (F_Tan, F_Normal)

#####################################################################################
#
# FUNCTION: scales the force matrices by desired percentiles of each in magnitude 
#           for colormap scalings 
#
#####################################################################################

def give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal):

    # Allocate Storage
    MagTan = np.zeros( (Npts,1) )
    MagNormal = np.zeros( (Npts,1) )

    # Find the magnitude of the force in each 
    for i in range(0,Npts):
        MagTan[i,0] = np.sqrt( F_Tan[i,0]*F_Tan[i,0]  + F_Tan[i,1]*F_Tan[i,1] )
        MagNormal[i,0]= np.sqrt( F_Normal[i,0]*F_Normal[i,0] + F_Normal[i,1]*F_Normal[i,1] )

    return (MagTan, MagNormal)

# # Finds Percentiles for Forces for Scaling (NOTE: This is in MATLAB syntax)
# prc90_T = prctile(MagTan,90);
# prc90_N = prctile(MagNormal,90);
# prc10_T = prctile(MagTan,10);
# prc10_N = prctile(MagNormal,10);
# 
# # "Scales the forces" via if-elseif statements by desired percentiles. 
# for i=1:Npts
#     
#     mT = MagTan(i);
#     mN = MagNormal(i);
#     
#     if mT >= prc90_T
#         MagTan(i) = prc10_T;
#     elseif mT <= prc10_T
#         MagTan(i) = prc10_T;
#     end
#     
#     if mN >= prc90_N
#         MagNormal(i) = prc10_N;
#     elseif mN <= prc10_N
#         MagNormal(i) = prc10_N;
#     end
#     
# end
