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
 
###############################################################################
#
# FUNCTION: creates the HEART-TUBE-EXAMPLE geometry 
#       and prints associated input files
#
################################################################################

from math import cos, sin, pi, sqrt
import numpy as np
import matplotlib.pyplot as plt

def HeartTube():

    #
    # Grid Parameters (MAKE SURE MATCHES IN input2d !!!)
    #
    Nx =  128        # # of Eulerian Grid Pts. in x-Direction (MUST BE EVEN!!!)
    Ny =  128        # # of Eulerian Grid Pts. in y-Direction (MUST BE EVEN!!!)
    Lx = 5.0        # Length of Eulerian Grid in x-Direction
    Ly = 5.0        # Length of Eulerian Grid in y-Direction


    # Immersed Structure Geometric / Dynamic Parameters #
    ds = Lx/(2*Nx)  # Lagrangian Spacing
    d = 1.0         # Diameter of the tube
    L = 3.0         # Length of the heart-tube
    struct_name = 'HeartTube' # Name for .vertex, .spring, etc files.


    # Call function to construct geometry
    xLag,yLag = give_Me_Immsersed_Boundary_Geometry(ds,L,d,Lx,Ly)


    # Plot Geometry to test
    plt.plot(xLag[:len(xLag)//2],yLag[:len(xLag)//2],'r-')
    plt.hold(True)
    plt.plot(xLag[len(xLag)//2:],yLag[len(xLag)//2:],'r-')
    plt.plot(xLag,yLag,'*')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis([0, Lx, 0, Ly])
    plt.show(block=True)


    # Prints .vertex file!
    print_Lagrangian_Vertices(xLag,yLag,struct_name)


    # Prints .spring file!
    k_Spring = 1e7
    print_Lagrangian_Springs(xLag,k_Spring,ds,struct_name)

    # Prints .muscle file! 
    #[a_f * Fmax *exp( -( (Q-1)/SK )^2 ) * (1/P0)*(b*P0-a*v)/(v+b); Q = LF/LFO]
    LFO = d; SK = 0.3; a = 0.25; b = 4.0; Fmax = 1e5
    print_Lagrangian_Muscles(xLag,LFO,SK,a,b,Fmax,struct_name)


    # Prints .beam file!
    k_Beam = 7.5e7; C = 0.0
    print_Lagrangian_Beams(xLag,k_Beam,C,struct_name)


    # Prints .target file!
    k_Target = 1e6
    print_Lagrangian_Target_Pts(xLag,k_Target,struct_name)

########################################################################
#
# FUNCTION: prints VERTEX points to a file called "struct_name".vertex
#
########################################################################

def print_Lagrangian_Vertices(xLag,yLag,struct_name):

    N = len(xLag)
    
    with open(struct_name + '.vertex','w') as vertex_fid:
        vertex_fid.write('{0}\n'.format(N))

        #Loops over all Lagrangian Pts.
        for s in range(N):
            vertex_fid.write('{0:1.16e} {1:1.16e}\n'.format(xLag[s],yLag[s]))

    
########################################################################
#
# FUNCTION: prints TARGET points to a file called "struct_name".target
#
########################################################################

def print_Lagrangian_Target_Pts(xLag,k_Target,struct_name):

    N = len(xLag)    # Total Number of Lagrangian Pts
    num = 1          # Number of target points on each end
    
    with open(struct_name + '.target','w') as target_fid:
        target_fid.write('{0}\n'.format(4*num))

        #Loops over all Lagrangian Pts.
        #for s in range(N)  
    
        #Left Bottom Target Points
        for s in range(num):
            target_fid.write('{0} {1:1.16e}\n'.format(s,k_Target))
        
        #Right Bottom Target Points
        for s in range(N//2-num,N//2):
            target_fid.write('{0} {1:1.16e}\n'.format(s,k_Target))
    
        #Left Top Target Points
        for s in range(N//2,N//2+1+(num-1)):
            target_fid.write('{0} {1:1.16e}\n'.format(s,k_Target))
    
        #Right Top Target Points
        for s in range(N-num,N):
            target_fid.write('{0} {1:1.16e}\n'.format(s,k_Target))

    
    
#####################################################################################
#
# FUNCTION: prints BEAM (Torsional Spring) points to a 
#        file called "struct_name".beam
#
#####################################################################################

def print_Lagrangian_Beams(xLag,k_Beam,C,struct_name):

    # k_Beam: beam stiffness
    # C: beam curvature
    
    N = len(xLag) # NOTE: Total number of beams = Number of Total Lag Pts. - 2
    
    with open(struct_name + '.beam','w') as beam_fid:
        beam_fid.write('{0}\n'.format(N-4))

        #spring_force = kappa_spring*ds/(ds**2)

        #BEAMS BETWEEN VERTICES
        for s in range(1,N-1):
            if s <= N//2 - 2:
                # Bottom of tube
                beam_fid.write('{0} {1} {2} {3:1.16e} {4:1.16e}\n'.format(
                    s-1,s,s+1,k_Beam,C))
            elif s >= N//2 + 1 and s <= N-2:
                # Top of tube
                beam_fid.write('{0} {1} {2} {3:1.16e} {4:1.16e}\n'.format(
                    s-1,s,s+1,k_Beam,C))
    

    
########################################################################
#
# FUNCTION: prints MUSCLE points to a file called "struct_name".muscle
#
########################################################################

def print_Lagrangian_Muscles(xLag,LFO,SK,a,b,Fmax,struct_name):

    N = len(xLag) #Number of Lagrangian Pts. Total
    
    with open(struct_name + '.muscle','w') as muscle_fid:
        muscle_fid.write('{0}\n'.format(N//2-2))

        #spring_force = kappa_spring*ds/(ds**2)

        #MUSCLES BETWEEN VERTICES
        for s in range(1,N//2-1):
            muscle_fid.write(
            '{0} {1} {2:1.16e} {3:1.16e} {4:1.16e} {5:1.16e} {6:1.16e}\n'.format(
            s, s+N//2, LFO, SK, a, b, Fmax))

    


########################################################################
#
# FUNCTION: prints SPRING points to a file called rubberband.spring
#
########################################################################

def print_Lagrangian_Springs(xLag,k_Spring,ds_Rest,struct_name):

    N = len(xLag) #Number of Lagrangian Pts. Total
    
    with open(struct_name + '.spring','w') as spring_fid:
        spring_fid.write('{0}\n'.format(N-2))

        #spring_force = kappa_spring*ds/(ds**2)

        #SPRINGS BETWEEN VERTICES
        for s in range(N):
            if s < N//2 - 1:
                spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(
                    s,s+1,k_Spring,ds_Rest))
            elif s >= N//2 and s < N-1:
                spring_fid.write('{0} {1} {2:1.16e} {3:1.16e}\n'.format(
                    s,s+1,k_Spring,ds_Rest))
    
    

########################################################################
#
# FUNCTION: creates the Lagrangian structure geometry
#
########################################################################

def give_Me_Immsersed_Boundary_Geometry(ds,L,d,Lx,Ly):
    
    # The immsersed structure is a straight heart-tube #
    import pdb; pdb.set_trace()
    x1 = list(np.arange(-L/2,L/2,ds))+[L/2] # Constructs x-Values for bottom of tube
    y1 = -d/2*np.ones(len(x1)) # Constructs y-Values for bottom of tube
    x1 = np.array(x1)

    x2 = np.array(x1)             # Constructs x-Values for top of tube
    y2 = -y1                      # Constructs y-Values for top of tube

    x1 = x1 + Lx/2                # Shift into correct box
    x2 = x2 + Lx/2                # Shift into correct box

    y1 = y1 + Ly/2                # Shift into correct box
    y2 = y2 + Ly/2                # Shift into correct box

    xLag = np.concatenate((x1, x2))
    yLag = np.concatenate((y1, y2))

    return (xLag,yLag)
    
if __name__ == "__main__":
    HeartTube()