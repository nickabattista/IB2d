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

function springs_info = update_Springs(dt,current_time,xLag,springs_info)

%springs_info: col 1: starting spring pt (by lag. discretization)
%              col 2: ending spring pt. (by lag. discretization)
%              col 3: spring stiffness
%              col 4: spring resting lengths


t = current_time;   %CURRENT TIME IN SIM.
ds1 = 3.906250e-03; %PHASE 1 Lag. Spacing
ds2 = 3.342919e-03; %PHASE 2 Lag. Spacing

Nbody = 106;                    %Number of springs before bell connections
RL = give_Resting_Lengths();    %Resting Lengths for Phase 1 (col 1) and Phase 2 (col 2)

% CHANGE RESTING LENGTH BTWN SIDES OF JELLYFISH BELL
for i=1:length( springs_info(:,4) )
    
    if i <= Nbody
        springs_info(i,4) = ds1 + 5*t*(ds2-ds1);
    else
        ii = i-Nbody;
        springs_info(i,4) = RL(ii,1) + 5*t*( RL(ii,2) - RL(ii,1) );
        %springs_info(107+(i-1),4) = distVec(i)*abs( cos(4*pi*current_time) );
        %springs_info(107+(i-1),4) = distVec(i) * ( 1 - 0.75*( sin(4*pi*current_time) ) );
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give Resting Lengths (hardcoded!)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RL = give_Resting_Lengths()

RL = [0.150000000000000   0.075000000000000
   0.149926741000240   0.074973170433545
   0.149706855857594   0.074892624017213
   0.149340020675537   0.074758187724185
   0.148825697410294   0.074569727429823
   0.148163136564709   0.074326763924451
   0.147351380996092   0.074028781527991
   0.146388643571571   0.073675421427274
   0.145274017365354   0.073265752667258
   0.144005152820154   0.072798927606369
   0.142580953464829   0.072274183901744
   0.140998544348621   0.071690309088121
   0.139255611062067   0.071045912762170
   0.137348541189626   0.070339929225114
   0.135274492318472   0.069570469632096
   0.133030516962324   0.068735403317456
   0.130612168717920   0.067832977122889
   0.128014600466564   0.066860709378195
   0.125234180745507   0.065816586316769
   0.122265614832230   0.064697078425658
   0.119103315657513   0.063499067862592
   0.115741451424401   0.062218444447808
   0.112173999448214   0.060850772863983
   0.108392818785359   0.059391311235482
   0.104393528844001   0.057835033327000
   0.100165653754868   0.056174749467196
   0.095702659718569   0.054402747848968
   0.090998237021989   0.052510837712120
   0.086041701055792   0.050487212358191
   0.080827140474999   0.048320369497431
   0.075349280754860   0.045994841644143
   0.069598549422225   0.043493120994599
   0.063573706279126   0.040793235757662
   0.057274941348722   0.037865974506005
   0.050704050274530   0.034676918247010
   0.043870119404512   0.031179908592683
   0.036787309077097   0.027316537892343
   0.029477909226729   0.023009349350403
   0.021975429586313   0.018167707789027
   0.014318966133393   0.012694598507574
   0.006558552552061   0.006558010776142];