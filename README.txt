README

-To Run Examples:
    1. Go into "Examples" Directory
    2. Choose which example you want to run and enter directory
    3. Type "main2d"
    4. You can change input data in the input2d data folder, or modify the geometry in the specified geometry file

-THIS CAN RUN 4 EXAMPLES:
    1. "Standard rubberband" - only uses springs between Lagrangian pts.
    2. "Wobbly Beam" - torsional springs between Lagrangian pts, w/ fixed ends via target points
    3. "Moving rubberband" - rubberband moving to the right, via updating target point positions
    4. "HeartTube" 
            ex1: - heart tube that pumps via muscle activation using combined length/tension - Hill model
            ex2: - heart tube that pumps via dynamic suction pumping by changing resting lengths of springs
            ex3: - heart tube that pumps via peristaltic waves of contraction by changing resting lengths of springs
    5. "Channel"
	    various channel flow examples with parabolic flow being produced by arbitrary force on Eulerian grid

-It has the ability to read in Lagrangian Point Data (.vertex), Springs (.spring), Torsional Springs (.beam),  
    Target Pts (.target), and Muscle Pts (.muscle)! 

-It has the function for updating target point positions, similar files for updating_springs() and update_beams() 
    can be made analogously. No file for update_beams() exists, yet.

-It can have pseudo-inflow conditions by inducing an arbitrary force onto the Eulerian grid

-It plots the Lagrangian structure along with Eulerian velocity, just the Lagrangian structure itself, 
    and a colormap of vorticity w/ immersed structure overlaid.

-NOTE: This code BLOWS UP when Lagrangian points cross a boundary!!!!!
