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
from math import cos, sin, pi, sqrt
import matplotlib.pyplot as plt
 
################################################################################
#
# FUNCTION: creates the RUBBERBAND-EXAMPLE geometry and prints associated input files
#
################################################################################

def Rubberband():

    #
    # Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
    #
    Nx =  64        # # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
    Ny =  64        # # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
    Lx = 1.0        # Length of Eulerian Grid in x-Direction
    Ly = 1.0        # Length of Eulerian Grid in y-Direction


    # Immersed Structure Geometric / Dynamic Parameters #
    N = 2*Nx        # Number of Lagrangian Pts. (2x resolution of Eulerian grid)
    ds = Lx/(2*Nx)  # Lagrangian Spacing
    a = 0.4         # Length of semi-major axis.
    b = 0.2         # Length of semi-minor axis.
    struct_name = 'rubberband' # Name for .vertex, .spring, etc files.


    # Call function to construct geometry
    xLag,yLag,C = give_Me_Immsersed_Boundary_Geometry(ds,N,a,b)


    # Plot Geometry to test
    plt.plot(xLag,yLag,'r-',xLag,yLag,'*')
    plt.xlabel('x'); plt.ylabel('y')
    plt.axis('square')
    plt.show(block=True)


    # Prints .vertex file!
    print_Lagrangian_Vertices(xLag,yLag,struct_name)


    # Prints .spring file!
    #k_Spring = 1e7
    #print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)


    # Prints .beam file!
    k_Beam = 1e7
    print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name)


    # Prints .target file!
    #k_Target = 1e7
    #print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

########################################################################
#
# FUNCTION: prints VERTEX points to a file called rubberband.vertex
#
########################################################################

def print_Lagrangian_Vertices(xLag,yLag,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.vertex', 'w') as vertex_fid:
        vertex_fid.write('{0}\n'.format(N))

        #Loops over all Lagrangian Pts.
        for s in range(N):
            X_v = xLag[s]
            Y_v = yLag[s]
            vertex_fid.write('{0:1.16e} {1:1.16e}\n'.format(X_v, Y_v))

    
########################################################################
#
# FUNCTION: prints Vertex points to a file called rubberband.vertex
#
########################################################################

def print_Lagrangian_Target_Pts(xLag,k_Target,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.target', 'w') as target_fid:
        target_fid.write('{0}\n'.format(N))

        #Loops over all Lagrangian Pts.
        for s in range(N):
            target_fid.write('{0:d} {1:1.16e}\n'.format(s, k_Target))

    
    
################################################################################
#
# FUNCTION: prints BEAM (Torsional Spring) points to a file called rubberband.beam
#
################################################################################

def print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name):

    # k_Beam: beam stiffness
    # C: beam curvature
    
    N = len(xLag) # NOTE: Total number of beams = Number of Total Lag Pts. - 2

    with open(struct_name + '.beam','w') as beam_fid:

        #v = range(0,N-1,2)
        
        beam_fid.write('{0:d}\n'.format(N))
        #beam_fid.write('{0:d}\n'.format(len(v)))

        #spring_force = kappa_spring*ds/(ds**2)

        #BEAMS BETWEEN VERTICES
        #for s = 1:2:N-1
        for s in range(N):
            if s==0:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                N-1, s, s+1, k_Beam, C[s] ))
            elif s <= N-2:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, s+1, k_Beam, C[s] ))
            else:
                #Case s=N-1
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, 0, k_Beam, C[s] ))


########################################################################
#
# FUNCTION: prints SPRING points to a file called rubberband.spring
#
########################################################################

def print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.spring', 'w') as spring_fid:

        spring_fid.write('{0:d}\n'.format(N))

        #spring_force = kappa_spring*ds/(ds**2)

        #SPRINGS BETWEEN VERTICES
        for s in range(N):
            if s < N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, s+1, k_Spring, ds_Rest))
            else:
                #Case s=N-1
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, 0, k_Spring, ds_Rest))


########################################################################
#
# FUNCTION: gives actually ELLIPTICAL piece
#
########################################################################

def compute_ELLIPTIC_Branch(ds,rmin,rmax):


    #initiate
    t = [0]
    xN = rmin*cos(0);   x = [xN]
    yN = rmax*sin(0);   y = [yN]

    while t[-1] <= 2*pi-ds: 
        
        xP = x[-1]           # x-Prev
        yP = y[-1]           # y-Prev
        
        tN = t[-1]            # previous angle
        tF = tN + pi/20         # far guess
        tGuess = (tN + tF)/2    # guess
        
        xN1 = rmin*cos(tGuess)   # x-guess 
        yN1 = rmax*sin(tGuess)   # y-guess
        err = ( ds - sqrt( (xN1-xP)**2 + (yN1-yP)**2 ) )
        
        while abs(err) > 1e-6:
           
            if err > 0:
              
                tN = tGuess              # Update 'close' PT. [tN,tGuess,tF]
                tGuess = (tN+tF)/2       # New Guess
                xN1 = rmin*cos(tGuess)   # x-guess 
                yN1 = rmax*sin(tGuess)   # y-guess          
              
            elif err < 0:
               
                tF = tGuess             # Update FAR PT. [tN,tGuess,tF] 
                tGuess = (tF+tN)/2      # New Guess
                xN1 = rmin*cos(tGuess)   # x-guess 
                yN1 = rmax*sin(tGuess)   # y-guess  
              
            
            #compute error
            err = ( ds - sqrt( (xN1-xP)**2 + (yN1-yP)**2 ) )


        #save values
        x.append(xN1)
        y.append(yN1)
        
        #update
        t.append(tGuess)


    #x.append(rmin*cos(angEnd))  
    #y.append(rmax*sin(angEnd))
        
    return (x,y,t)
    
    
    
########################################################################
#
# FUNCTION: creates the Lagrangian structure geometry
#
########################################################################

def give_Me_Immsersed_Boundary_Geometry(ds,N,rmin,rmax):

    # The immsersed structure is an ellipse #
    xLag,yLag,angs = compute_ELLIPTIC_Branch(ds,rmin,rmax)
    xLag = [item+0.5 for item in xLag]
    yLag = [item+0.5 for item in yLag]

    # COMPUTES CURAVTURE IF WANT TO STAY IN INITIAL CONFIGURATION
    C = compute_Curvatures(ds,angs,rmin,rmax,xLag,yLag)

    N = len(xLag)
    r_eff = sqrt(rmin*rmax)
    xLag2 = []; yLag2 = []
    for ii in range(N):
        
        xLag2.append(0.5 + r_eff * cos( 2*pi/N*ii ))
        yLag2.append(0.5 + r_eff * sin( 2*pi/N*ii ))

    # COMPUTES CURAVTURE IF WANT TO SETTLE INTO A CIRCLE
    C = compute_Curvatures(ds,angs,rmin,rmax,xLag2,yLag2)

    return (xLag,yLag,C)


########################################################################
#
# FUNCTION: computes "curvature" of ellipse
# 
# NOTE: not curvature in the traditional geometric sense, in the 'discrete'
# sense through cross product.
#
########################################################################

def compute_Curvatures(ds,angs,rmin,rmax,xLag,yLag):

    #a-x component (rmin)
    #b-y component (rmax)
    #C = ab / ( sqrt( a^2*sin(t)^2 + b^2*cos(t)^2  )  )^3

    N = len(xLag)
    C = np.zeros(len(angs))

    #Note: -needs to be done same order as you print .beam file!
    #      -THIS MAKES INITIAL BEAM CONFIGURATION THE DESIRED CURAVTURE!!
    for ii in range(N):
       
        if ii==0:
       
            # Pts Xp -> Xq -> Xr (same as beam force calc.)
            Xp = xLag[-1]; Xq = xLag[ii]; Xr = xLag[ii+1]
            Yp = yLag[-1]; Yq = yLag[ii]; Yr = yLag[ii+1]
       
        elif ii<N-1:
       
            # Pts Xp -> Xq -> Xr (same as beam force calc.)
            Xp = xLag[ii-1]; Xq = xLag[ii]; Xr = xLag[ii+1]
            Yp = yLag[ii-1]; Yq = yLag[ii]; Yr = yLag[ii+1]
           
        else:
           
            # Pts Xp -> Xq -> Xr (same as beam force calc.)
            Xp = xLag[ii-1]; Xq = xLag[ii]; Xr = xLag[0]
            Yp = yLag[ii-1]; Yq = yLag[ii]; Yr = yLag[0]
        
        # Small numbers here, roundoff error can make result slightly different
        #   from what you get in MATLAB.
        C[ii] = (Xr-Xq)*(Yq-Yp) - (Yr-Yq)*(Xq-Xp) #Cross product btwn vectors
       
    return C

if __name__ == "__main__":
    Rubberband()