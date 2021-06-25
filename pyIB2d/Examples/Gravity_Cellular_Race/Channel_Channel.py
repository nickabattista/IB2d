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
# FUNCTION: creates the CHANNEL_CHANNEL-EXAMPLE geometry and prints associated input files
#
################################################################################

def Channel_Channel():

    #
    # Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
    #
    Nx =  64        # # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
    Ny =  64        # # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
    Lx = 1.0        # Length of Eulerian Grid in x-Direction
    Ly = 1.0        # Length of Eulerian Grid in y-Direction


    # Immersed Structure Geometric / Dynamic Parameters #
    ds = min(Lx/(2*Nx),Ly/(2*Ny))   # Lagrangian spacing
    L = 0.9*Lx                      # Length of Channel
    w = 0.2*Ly                      # Width of Channel
    x0 = Lx/2                       # x-Center for Cylinder
    y0 = 0.8                        # y-Center for Cylinder
    r = w/6                         # Radii of Cylinder
    struct_name = 'channel' # Name for .vertex, .spring, etc files.


    # Call function to construct geometry
    xLag_C,yLag_C = give_Me_Cylinder_Immersed_Boundary_Geometry(ds,r,x0,y0)
    
    #Left Cell
    xLag_C1 = xLag_C-0.25; yLag_C1 = yLag_C
    
    #Middle Cell
    xLag_C2 = xLag_C; yLag_C2 = yLag_C

    #Right Cell
    xLag_C3 = xLag_C+0.25; yLag_C3 = yLag_C

    # Plot Geometry to test
    plt.plot(xLag_C1,yLag_C1,'r-')
    plt.plot(xLag_C2,yLag_C2,'r-')
    plt.plot(xLag_C3,yLag_C3,'r-')
    #
    plt.plot(xLag_C1,yLag_C1,'g*')
    plt.plot(xLag_C2,yLag_C2,'g*')
    plt.plot(xLag_C3,yLag_C3,'g*')
    #
    plt.xlabel('x'); plt.ylabel('y')
    plt.axis([0, Lx, 0, Ly])
    plt.show(block=True)
    
    # Combine all lagrangian pts into one vector
    # Add xLagPts from Circle to xLag Pt. Vector (*no springs or beams*)
    xLag = np.hstack((xLag_C1, xLag_C2, xLag_C3))
    # Add yLagPts from Circle to yLag Pt. Vector (*no springs or beams*)
    yLag = np.hstack((yLag_C1, yLag_C2, yLag_C3)) 
    
    # Prints .vertex file!
    print_Lagrangian_Vertices(xLag,yLag,struct_name)

    # Prings .mass file!
    k_Mass = 1e6                # 'spring' stiffness parameter for tethering
    Mass = [5e-2, 2e-1, 1e0]    # "MASS" value for 'ghost' nodal movement
    print_Lagrangian_Mass_Pts(xLag_C1,k_Mass,Mass,struct_name)
    
    # Prints .spring file!
    k_Spring = 5e7
    ds_Rest = sqrt( (xLag_C[0]-xLag_C[1])**2 + (yLag_C[0]-yLag[1])**2 )
    print_Lagrangian_Springs(xLag_C,yLag_C,k_Spring,ds_Rest,r,struct_name)


    # Prints .beam file!
    k_Beam = 1e7; C = 0.0 #1/r      # curvature of a circle
    print_Lagrangian_Beams(xLag_C,yLag_C,k_Beam,C,struct_name)


    # Prints .target file!
    # k_Target = 1e7; Noff = 3*xLag_C.size
    # print_Lagrangian_Target_Pts(xLag_T,Noff,k_Target,struct_name)

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
# FUNCTION: prints TARGET points to a file called rubberband.target
#
########################################################################

def print_Lagrangian_Target_Pts(xLag,Noff,k_Target,struct_name):

    N = 3*xLag.size # for 3 channels
    
    with open(struct_name + '.target', 'w') as target_fid:
        target_fid.write('{0}\n'.format(N))

        #Loops over all Lagrangian Pts.
        for s in range(N):
            # choose correct lag. pt.
            target_fid.write('{0:d} {1:1.16e}\n'.format(Noff+s, k_Target))

            
########################################################################
#
# FUNCTION: prints MASS points to a file called struct_name.mass
#
########################################################################  
    
def print_Lagrangian_Mass_Pts(xLag,kMass,Mass,struct_name):

    #LOOP OVER LAG PTS FOR MASS PTS IN LAGRANGIAN INDEXING

    N = xLag.size
    
    with open(struct_name + '.mass','w') as mass_fid:
        mass_fid.write('{0}\n'.format(3*N))

        #Loops over all Lagrangian Pts.
        for s in range(N):
            mass_fid.write('{0:d} {1:1.16e} {2:1.16e}\n'.format(s,kMass,Mass[0]))
        for s in range(N,2*N):
            mass_fid.write('{0:d} {1:1.16e} {2:1.16e}\n'.format(s,kMass,Mass[1]))
        for s in range(2*N,3*N):
            mass_fid.write('{0:d} {1:1.16e} {2:1.16e}\n'.format(s,kMass,Mass[2]))

        
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

        beam_fid.write('{0:d}\n'.format(3*N))

        #spring_force = kappa_spring*ds/(ds**2)

        #BEAMS BETWEEN VERTICES
        for s in range(N):
            if s==0:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                N-1, s, s+1, k_Beam, C))
            elif s <= N-2:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, s+1, k_Beam, C))
            else:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, 0, k_Beam, C))
        
        for s in range(N,2*N):
            if s==N:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                2*N-1, s, s+1, k_Beam, C))
            elif s <= 2*N-2:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, s+1, k_Beam, C))
            else:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, N, k_Beam, C))
                
        for s in range(2*N,3*N):
            if s==2*N:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                3*N-1, s, s+1, k_Beam, C))
            elif s <= 3*N-2:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, s+1, k_Beam, C))
            else:
                beam_fid.write('{0:d} {1:d} {2:d} {3:1.16e} {4:1.16e}\n'.format(\
                s-1, s, 2*N, k_Beam, C))


########################################################################
#
# FUNCTION: prints SPRING points to a file called rubberband.spring
#
########################################################################

def print_Lagrangian_Springs(xLag,yLag,k_Spring,ds_Rest,r,struct_name):

    N = xLag.size
    
    with open(struct_name + '.spring', 'w') as spring_fid:

        spring_fid.write('{0}\n'.format(3*N + 3/2*N)) #N MUST BE EVEN!

        #spring_force = kappa_spring*ds/(ds**2)

        #SPRINGS BETWEEN VERTICES
        for s in range(3*N):
            if s < N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, s+1, k_Spring, ds_Rest))
            elif s==N-1:
                #Case s=N-1
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, 0, k_Spring, ds_Rest))
            elif s >= N and s < 2*N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, s+1, k_Spring, ds_Rest))
            elif s==2*N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, N, k_Spring, ds_Rest))
            elif s >= 2*N and s <3*N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, s+1, k_Spring, ds_Rest))
            elif s==3*N-1:
                spring_fid.write('{0:d} {1:d} {2:1.16e} {3:1.16e}\n'.format(\
                s, 2*N, k_Spring, ds_Rest))
                
        for s in range(int(N/2)):
            spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(\
                s, int(s+N/2), k_Spring/100, 2*r))
        for s in range(int(N/2)):
            spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(\
                N+s, int(N+s+N/2), k_Spring/100, 2*r))
        for s in range(int(N/2)):
            spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(\
                2*N+s, int(2*N+s+N/2), k_Spring/100, 2*r))

                
########################################################################
#
# FUNCTION: creates the Lagrangian structure geometry for channel
#
########################################################################

def give_Me_Channel_Immsersed_Boundary_Geometry(ds,L,w,Lx,Ly):

    # The immsersed structure is a channel #
    x = np.arange(Lx-L,L+(Lx-L)/2,ds)   #xPts
    yBot = (Ly-w)/2 #yVal for bottom of Channel
    yTop = Ly - (Ly-w)/2 #yVal for top of Channel
    
    xLag = np.concatenate((x,x))
    yLag = np.concatenate((yBot*np.ones(x.size), yTop*np.ones(x.size)))
    
    return (xLag,yLag)
    
    
########################################################################
#
# FUNCTION: creates the Lagrangian structure geometry for cylider
#
########################################################################

def give_Me_Cylinder_Immersed_Boundary_Geometry(ds,r,x0,y0):
    
    # The immersed structure is a cylinder #
    
    dtheta = ds/(2*r)
    theta = 0
    xLag = []; yLag = []
    while theta < 2*pi:
        xLag.append(x0 - r*cos(theta))
        yLag.append(y0 - r*sin(theta))
        theta += dtheta
    xLag = np.array(xLag)
    yLag = np.array(yLag)
    
    return (xLag,yLag)

if __name__ == "__main__":
    Channel_Channel()