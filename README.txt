IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
 	fluid-structure interaction models. This version of the code is based off of
	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.

Author: Nicholas A. Battista
Email:  nick.battista@unc.edu
Date Created: May 27th, 2015
Institution: UNC-CH

This code is capable of creating Lagrangian Structures using:
 	1. Springs
 	2. Beams (*torsional springs)
 	3. Target Points
	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")

One is able to update those Lagrangian Structure Parameters, e.g., spring constants, resting lengths, etc
 
There are a number of built in Examples, mostly used for teaching purposes. 

If you would like us to add a specific muscle model, please contact Nick (nick.battista@unc.edu) 

If you use this code for the purposes of teaching, please let Nick know as well :)

%*****************************************************************************%

%*****************************% HAPPY COMPUTING! %****************************%

%*****************************************************************************%

-To Run Examples:
    1. Go into "Examples" Directory
    2. Choose which example you want to run and enter directory
    3. Type "main2d"
    4. You can change input data in the input2d data folder, or modify the geometry in 
       the specified geometry file

-THIS CAN RUN 5 DIFFERENT EXAMPLES:
    1. "Standard Rubberband" - only uses springs between Lagrangian pts.
    2. "Wobbly Beam" - torsional springs between Lagrangian pts, w/ fixed ends 
        via target points
    3. "Moving Rubberband" - rubberband moving to the right, via updating target 
        point positions
    4. "HeartTube" 
            ex1: - heart tube that pumps via muscle activation using combined 
                   length/tension - Hill model
            ex2: - heart tube that pumps via dynamic suction pumping by changing 
                   resting lengths of springs
            ex3: - heart tube that pumps via peristaltic waves of contraction by 
                   changing resting lengths of springs
    5. "Channel_Flow" - various channel flow examples with parabolic flow being 
                        produced by arbitrary force on Eulerian grid
    6. "Pulsing_Heart" - cartoon heart that "pumps" via moving target points 
                         (interpolates between two phases)

-It has the ability to read in Lagrangian Point Data (.vertex), Springs (.spring), 
    Torsional Springs (.beam),  Target Pts (.target), and Muscle Pts (.muscle)! 

-It has the capabilitiy for updating model data:
    a. updating target point positions -> Example in Moving_Rubberband & Pulsing_Heart
    b. updating_springs() -> Examples in HeartTube
    c. update_beams() [can be made analogously] 
    d. update_muscles() [can be made analogously]

-It can have pseudo-inflow conditions by inducing an arbitrary force onto the Eulerian 
    grid (e.g., Channel Flow Examples)

-It can plot the following things in Matlab, if plot_Matlab flag = 1 in input2d:
    a. Vorticity (colormap) + Lagrangian Pts.
    b. Magnitude Velocity (colormap) + Lagrangian Pts.
    c. Pressure (colormap) + Lagrangian Pts.
    d. Velocity (vector form) + Lagrangian Pts.
    e. Lagrangian Pts. themselves 

-It has a flag for print dump interval (shared between printing to .vtk format + Matlab 
    plotting)

-NOTE: This code BLOWS UP when Lagrangian points cross a boundary!!!!!
