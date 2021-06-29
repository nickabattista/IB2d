NOTE:

[1] Before running main2d.m, please run Make_Tub_Geometry_and_Initial_Concentration.m
        to create the input files for the simulation! :)

[2] This simulation has 4 (horizontal) regions of horizontal flow (every other moving in same direction)
        to create shear at the interfaces. There is then vertical perturbations in the flow
        as prescribed in "please_Compute_External_Forcing.m" to create instability

[2] This example doesn't use the Boussinesq approximation - just for illustrative purposes 
        of what the KH-Instability qualitatively looks like

[3] Simulation at 512x1024 resolution takes about 3.5 hours to run on 2.7 GHz processor
