%*****************************************************************************%
%*******************************% IB2d Python %*******************************%
%*****************************************************************************%

IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based 
	off of Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

Author: Nicholas A. Battista
Email:  nick.battista@unc.edu
Date Created: May 27th, 2015
Institution: University of North Carolina at Chapel Hill
Website: http://battista.web.unc.edu
GitHub: http://www.github.com/nickabattista

Python 3.5 port by: Christopher Strickland
Date Created: March 3rd, 2016
Institution: University of North Carolina at Chapel Hill
Website: http://www.christopherstrickland.info
GitHub: http://www.github.com/mountaindust

Dependencies:
 	- NumPy
 	- Matplotlib
 	- Numba

Additionally, an optional C library (plus Cython glue code) is provided to improve the speed at which data writes out to disk. For info on compiling this library, please see Cython_README.md in the IBM_Blackbox directory.

This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
	4. Muscle-Model 1 (combined Force-Length-Velocity model with a 
            Hill (i.e.,Length-Tension) Model )
    5. 3-Element Hill Muscle (combined 3-element hill mode coupled w/
            Force-Velocity/Length-Tension Model) -- UNTESTED!! THIS IS LIKELY BROKEN.
    6. Mass Points (with or without influence of gravity)
    7. Porous Structures

One is able to update those Lagrangian Structure Parameters, e.g., spring constants, resting lengths, etc.
 
There are a number of ported Examples, mostly used for testing and teaching purposes. Many more are available in the Matlab code - if you port another example, please send us a pull request!

If you would like us to add a specific muscle model, please contact Nick (nick.battista@unc.edu) 

If you use this code for the purposes of teaching, research, or recreation please let Nick know as well :)

%*****************************************************************************%
%*****************************% HAPPY COMPUTING! %****************************%
%*****************************************************************************%

-To Run Examples:
    1. Go into "Examples" Directory
    2. Choose which example you want to run and enter directory
    3. Type "python main2d.py"
    4. You can change input data in the input2d data folder, or modify the 
       geometry in the specified geometry file

    Alternatively you can run pyIB2d.py, specifying a directory containing
    the simulation files via one of its flags. Type "python pyIB2d.py -h" for
    further info. Note: with this method, all data will be written into this
    folder (IB2d/pyIB2d) rather than the simulation folder.

-THIS VERSION HAS 6 DIFFERENT FLAVORS OF EXAMPLES:
    1. "Standard Rubberband" - only uses springs/beams between Lagrangian pts.
    2. "HeartTube" - example of pumping in tubular hearts	
    3. "Tracers" - example of inserting tracer particles into simulations
            Tracers in an impedance pump
    4. "Mass_Points" - "adding artificial mass" to the Lagrangian structure
            'massive' cells racing due to underlying fluid velocity only!
    5. "Porous_Rubberband" - added porosity to standard rubberband problem
    6. "Concentration_Diffusion" - Pure Diffusion of Passive Scalar in Box

-It has the ability to read in Lagrangian Point Data (.vertex), Springs 
	(.spring), Torsional Springs (.beam),  Target Pts (.target), and 
	Muscle Pts (.muscle), 3-Hill-Muscle-Pts (.muscle_Hill), 
        Tracer Particles (.tracer), Mass Points (.mass), and  Porous Media 
        (.porous), and initial concentration (.concentration)! 

-It has the capabilitiy for updating model data

-It can have pseudo-inflow conditions by inducing an arbitrary force onto the 
    Eulerian grid (e.g., Channel Flow Examples)

-You can choose to have gravity exerting forced (yes/no) as well as prescribe
    the direction of gravitational influence in any direction you'd wish, 
    in input2d

-You can have a background concentration gradient that is advected and diffused
    via the background flow.

-It can plot the following things in matplotlib, if plot_Matlab flag = 1 in input2d:
    a. Vorticity (colormap) + Lagrangian Pts.
    b. Magnitude Velocity (colormap) + Lagrangian Pts.
    c. Pressure (colormap) + Lagrangian Pts.
    d. Velocity (vector form) + Lagrangian Pts.
    e. Lagrangian Pts. themselves 

-It has a flag for print dump interval (shared between printing to .vtk format 
    + matplotlib plotting)

-NOTE: This code *may* BLOW UP when Lagrangian points cross a boundary!!!!!

%*****************************************************************************%
%*******************************% VISUALIZATION %*****************************%
%*****************************************************************************%

-These examples print data as .vtk files, which can be read by Paraview and 
    VisIt.

-Every example prints the following:

        LAGRANGIAN PTS:     a. Lag. Pts. themselves
                            b. Lag. Pts. w/ connections!

        SCALARS (colormap): a. Vorticity
                            b. Magnitude of Velocity
                            c. uX (x-directed velocity)
                            d. uY (y-directed velocity)
                            e. Pressure

        VECTORS: a. velocity data

        TRACERS: a. tracer particle locations

-There are flags in input2d for the printing interval between saving data
        -> print_Dump

-There are flags in input2d whether you'd like Matlab to plot various 
    quantities as simulation progresses:

            -> plot_Matlab     (set = 1, if yes, have Matlab plot)
            -> plot_<Quantity> (set = 1, if yes for quantity)
            
