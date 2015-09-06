%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: updates the spring attributes!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function beams_info = update_Beams(dt,current_time,beams_info)


% beams_info:   col 1: 1ST PT.
%               col 2: MIDDLE PT. (where force is exerted)
%               col 3: 3RD PT.
%               col 4: beam stiffness
%               col 5: curavture


period = 0.1;                    % Period of contraction/expansion
P1 = 0.3*period;                 % Fraction of Period for Contraction
P2 = period - P1;                % Fraction of Period for Expansion

t = rem( current_time, period);  % Current time in simulation, via modular arithmetic w/ period

C = give_Beam_Curvatures();      % Gives beam CURVATURES for each phase


% CHANGE RESTING LENGTH BTWN SIDES OF JELLYFISH BELL
for i=1:length( beams_info(:,5) )
    
    if t<= P1
        beams_info(:,5) = C(:,1) + (t/P1)*( C(:,2) - C(:,1) );
    else
        tt = t - P1;
        beams_info(:,5) = C(:,2) + ( tt/P2 )*( C(:,1) - C(:,2) );
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give Beam (Torsional Spring) CURVATURES (hardcoded!)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function curves = give_Beam_Curvatures()

curves = [-2.864292e-07 -8.979307e-08
-2.872202e-07 -9.008286e-08
-2.885419e-07 -9.047027e-08
-2.903989e-07 -9.105332e-08
-2.927983e-07 -9.193269e-08
-2.960176e-07 -9.291650e-08
-2.998039e-07 -9.410602e-08
-3.041696e-07 -9.560733e-08
-3.091325e-07 -9.727721e-08
-3.147086e-07 -9.917227e-08
-3.212043e-07 -1.014054e-07
-3.286533e-07 -1.038846e-07
-3.370976e-07 -1.066737e-07
-3.462828e-07 -1.099535e-07
-3.565460e-07 -1.136386e-07
-3.682544e-07 -1.177511e-07
-3.811728e-07 -1.223174e-07
-3.953674e-07 -1.274862e-07
-4.112494e-07 -1.334263e-07
-4.289288e-07 -1.401956e-07
-4.485268e-07 -1.479964e-07
-4.705421e-07 -1.569195e-07
-4.947855e-07 -1.670684e-07
-5.218025e-07 -1.788625e-07
-5.522347e-07 -1.928000e-07
-5.855227e-07 -2.091336e-07
-6.227644e-07 -2.286900e-07
-6.643237e-07 -2.523007e-07
-7.096185e-07 -2.809690e-07
-7.599430e-07 -3.163939e-07
-8.152455e-07 -3.607256e-07
-8.742746e-07 -4.178003e-07
-9.372455e-07 -4.927781e-07
-1.003177e-06 -5.936574e-07
-1.070229e-06 -7.336370e-07
-1.136246e-06 -9.341778e-07
-1.197397e-06 -1.228301e-06
-1.250119e-06 -1.661816e-06
-1.290668e-06 -2.263332e-06
-1.016078e-06 -2.875063e-06
-7.840238e-07 -3.139832e-06
-1.016078e-06 -2.875063e-06
-1.290668e-06 -2.263332e-06
-1.250119e-06 -1.661816e-06
-1.197397e-06 -1.228301e-06
-1.136246e-06 -9.341778e-07
-1.070229e-06 -7.336370e-07
-1.003177e-06 -5.936574e-07
-9.372455e-07 -4.927781e-07
-8.742746e-07 -4.178003e-07
-8.152455e-07 -3.607256e-07
-7.599430e-07 -3.163939e-07
-7.096185e-07 -2.809690e-07
-6.643237e-07 -2.523007e-07
-6.227644e-07 -2.286900e-07
-5.855227e-07 -2.091336e-07
-5.522347e-07 -1.928000e-07
-5.218025e-07 -1.788625e-07
-4.947855e-07 -1.670684e-07
-4.705421e-07 -1.569195e-07
-4.485268e-07 -1.479964e-07
-4.289288e-07 -1.401956e-07
-4.112494e-07 -1.334263e-07
-3.953674e-07 -1.274862e-07
-3.811728e-07 -1.223174e-07
-3.682544e-07 -1.177511e-07
-3.565460e-07 -1.136386e-07
-3.462828e-07 -1.099535e-07
-3.370976e-07 -1.066737e-07
-3.286533e-07 -1.038846e-07
-3.212043e-07 -1.014054e-07
-3.147086e-07 -9.917227e-08
-3.091325e-07 -9.727721e-08
-3.041696e-07 -9.560733e-08
-2.998039e-07 -9.410602e-08
-2.960176e-07 -9.291650e-08
-2.927983e-07 -9.193269e-08
-2.903989e-07 -9.105332e-08
-2.885419e-07 -9.047027e-08
-2.872202e-07 -9.008286e-08
-2.864292e-07 -8.979307e-08];