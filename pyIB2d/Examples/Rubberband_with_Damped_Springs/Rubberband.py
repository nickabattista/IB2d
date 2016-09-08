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
    Nx =  32        # # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
    Ny =  32        # # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
    Lx = 1.0        # Length of Eulerian Grid in x-Direction
    Ly = 1.0        # Length of Eulerian Grid in y-Direction


    # Immersed Structure Geometric / Dynamic Parameters #
    N = 2*Nx         # Number of Lagrangian Pts. (2x resolution of Eulerian grid)
    a = 0.4          # Length of semi-major axis.
    b = 0.2          # Length of semi-minor axis.
    ds_Rest = 0      # Resting length of springs
    struct_name = 'rubberband' # Name for .vertex, .spring, etc files.


    # Call function to construct geometry
    xLag,yLag = give_Me_Immsersed_Boundary_Geometry(N,a,b)


    # Prints .vertex file!
    print_Lagrangian_Vertices(xLag,yLag,struct_name)


    # Prints .spring file!
    # k_Spring = 2.5e4
    # print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,struct_name)

    # Prints .d_spring file!
    k_Spring = 2.5e4
    b_damp = 5.0
    print_Lagrangian_Damped_Springs(xLag,yLag,k_Spring,ds_Rest,b_damp,struct_name)


    # Prints .beam file!
    #k_Beam = 0.5; C = 0.0
    #print_Lagrangian_Beams(xLag,yLag,k_Beam,C,struct_name)


    # Prints .target file!
    #k_Target = 1e7
    #print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)


    # Plot Geometry to test
    plt.plot(xLag,yLag,'r-',xLag,yLag,'*')
    plt.xlabel('x'); plt.ylabel('y')
    plt.axis('square')
    plt.show(block=True)


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

        beam_fid.write('{0:d}\n'.format(N))

        #spring_force = kappa_spring*ds/(ds**2)

        #BEAMS BETWEEN VERTICES
        for s in range(N):
            if  s <= N-1:       
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


##############################################################################
#
# FUNCTION: prints DAMPED SPRING points to a file called rubberband.d_spring
#
##############################################################################

def print_Lagrangian_Damped_Springs(xLag,yLag,k_Spring,ds_Rest,b_damp,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.d_spring', 'w') as spring_fid:

        spring_fid.write('{0:d}\n'.format(N))

        #spring_force = kappa_spring*ds/(ds**2)

        #SDAMPED PRINGS BETWEEN VERTICES
        for s in range(N):
            if s < N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e} {4:1.16e}\n'.format(\
                s, s+1, k_Spring, ds_Rest, b_damp))
            else:
                #Case s=N-1
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e} {4:1.16e}\n'.format(\
                s, 0, k_Spring, ds_Rest, b_damp))                

########################################################################
#
# FUNCTION: creates the Lagrangian structure geometry
#
########################################################################


def give_Me_Immsersed_Boundary_Geometry(N,rmin,rmax):

    # The immsersed structure is an ellipse #
    t = 2*pi/N
    xLag = np.zeros(N)
    yLag = np.zeros(N)
    for ii in range(N):
    
        xLag[ii] = 0.5 + rmax * cos( 2*pi/N*ii )
        yLag[ii] = 0.5 + rmin * sin( 2*pi/N*ii )
    
    return (xLag,yLag)

if __name__ == "__main__":
    Rubberband()