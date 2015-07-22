%*****************************************************************************%
%***********************************% IB2d %**********************************%
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

This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
	4. Muscle-Model (combined Force-Length-Velocity model, "Hill(Length-
                        Tension)")
        5. Mass Points (with or without influence of gravity)

One is able to update those Lagrangian Structure Parameters, e.g., spring constants, resting lengths, etc
 
There are a number of built in Examples, mostly used for teaching purposes. 

If you would like us to add a specific muscle model, please contact Nick (nick.battista@unc.edu) 

If you use this code for the purposes of teaching, research, or recreation please let Nick know as well :)

%*****************************************************************************%
%*****************************% HAPPY COMPUTING! %****************************%
%*****************************************************************************%

-To Run Examples:
    1. Go into "Examples" Directory
    2. Choose which example you want to run and enter directory
    3. Type "main2d"
    4. You can change input data in the input2d data folder, or modify the 
       geometry in the specified geometry file

-THIS VERSION HAS 8 DIFFERENT FLAVORS OF EXAMPLES:
    1. "Standard Rubberband" - only uses springs between Lagrangian pts.
    2. "Wobbly Beam" - torsional springs between Lagrangian pts, w/ fixed ends 
        via target points
    3. "Moving Rubberband" - rubberband moving along a zig-zag pattern, via 
        updating target point positions
    4. "HeartTube" - tubular heart pumping examples
            ex1: - heart tube that pumps via muscle activation using combined 
                   length/tension - Hill model
            ex2: - heart tube that pumps via dynamic suction pumping by 
                   changing resting lengths of springs
            ex3: - heart tube that pumps via peristaltic waves of contraction 
                   by changing resting lengths of springs
    5. "Lymphangion" - tube that pumps via muscle activation using combined a
                   hybrid length/tension and Hill model
    6. "Channel_Flow" - various channel flow examples with parabolic flow being 
                        produced by arbitrary force on Eulerian grid
            ex1: - Flow around a cylinder in a channel
            ex2: - Flow around a cylinder with only one wall of channel
            ex3: - Flow in a channel with a hole in the channel
            ex4: - Flow produced in the middle of the channel
            ex5: - Flow produced in an empty channel
            ex6: - Flow in a channel, mocking an aneurism
    7. "Pulsing_Heart" - cartoon heart that "pumps" via moving target points 
                         (interpolates between two phases)
    8. "Tracers" - examples of inserting tracer particles into simulations
            ex1: - Tracers in channel flow 
            ex2: - Tracers in an impedance pump
            ex3: - Tracers in a peristaltic pump
    9. "Mass_Points" - 'massive' cell in pulsatile channel flow
                     - 'massive' cells racing due to underlying fluid velocity!
                     - 'massive' cells racing under gravity
                     - 'massive' cells in 'gravity vs. pulsating flow'

-It has the ability to read in Lagrangian Point Data (.vertex), Springs 
	(.spring), Torsional Springs (.beam),  Target Pts (.target), and 
	Muscle Pts (.muscle), Tracer Particles (.tracer), Mass Points (.mass)! 

-It has the capabilitiy for updating model data:
    a. updating target point positions -> Example in Moving_Rubberband & 
       Pulsing_Heart
    b. updating_springs() -> Examples in HeartTube
    c. update_beams() [can be made analogously] 
    d. update_muscles() [can be made analogously]

-It can have pseudo-inflow conditions by inducing an arbitrary force onto the 
    Eulerian grid (e.g., Channel Flow Examples)

-You can choose to have gravity exerting forced (yes/no) as well as prescribe
    the direction of gravitational influence in any direction you'd wish, 
    in input2d

-It can plot the following things in Matlab, if plot_Matlab flag = 1 in input2d:
    a. Vorticity (colormap) + Lagrangian Pts.
    b. Magnitude Velocity (colormap) + Lagrangian Pts.
    c. Pressure (colormap) + Lagrangian Pts.
    d. Velocity (vector form) + Lagrangian Pts.
    e. Lagrangian Pts. themselves 

-It has a flag for print dump interval (shared between printing to .vtk format 
    + Matlab plotting)

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
            
