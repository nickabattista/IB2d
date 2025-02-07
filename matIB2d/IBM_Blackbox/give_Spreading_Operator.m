%------------------------------------------------------------------------------------------
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled non-linear 
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn[@]tcnj[.]edu
% 
% IB2d was Created: May 27th, 2015 at UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
% 	4. Muscle-Model (combined Force-Length-Velocity model, "Hill+(Length-Tension)")
%   .
%   .
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
%
%------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Returns the SPREADING operator for spreading information
%           from the Lagrangian Grid --> Eulerian Grid
% 
%      INPUTS:
%           [1] xLag,yLag: Lagrangian Pt Values in x and y
%           [2]     dx,dy: Grid resolution in x and y
%           [3]     NLags: # of Lagrangian Points
%           [4]     Nx,Ny: # of grid cells in x and y
%
%      OUTPUT: S: Spreading Operator (of size S\n R^{(Nx*Ny) x N} )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = give_Spreading_Operator(xLag,yLag,dx,dy,NLags,Nx,Ny)
 
    %----------------------------------------------------------------
    % Allocate storage for the Spreading Operator, S
    %    --> of Size (Nx*Ny) x NLags
    %    --> Each column will correspond to a different Lag. Pt. 
    %    --> Each column represents all points in Eulerian Grid
    %        (it gets flattened from a 2-D array to 1-D column vector)
    %    --> Since using 4-pt stencil, only need 16 non-zero entries
    %        in any given column
    %    --> Therefore, S has 16*NLags possible non-zero entries
    %    --> Spalloc(m,n,nz) creates an all-zero sparse matrix of  
    %        size m-by-n with room to hold nz nonzero elements.
    %----------------------------------------------------------------
    S = spalloc( Nx*Ny, NLags, 16*NLags);

    %-------------------------------------------  
    % Loop over each Lagrangian point
    %-------------------------------------------
    for n=1:NLags

        %------------------------------------------------------   
        % Scale n-th Lag. Point to help find Eulerian Grid ID
        %------------------------------------------------------ 
        xAux = xLag(n)/dx + 1;
        yAux = yLag(n)/dy + 1;

        %------------------------------------------------------
        % Find the next lowest Eulerian grid ID
        %   (so Lag Pt. should be "up and to the right")
        %------------------------------------------------------
        idX_Seed = floor( xAux );
        idY_Seed = floor( yAux );

        %------------------------------------------------------
        % Get Eulerian grid indices to spread across x and y
        %       (assumes a 4-pt delta function stencil)
        %------------------------------------------------------
        indsX = idX_Seed-1:idX_Seed+2;
        indsY = idY_Seed-1:idY_Seed+2;

        %-----------------------------------------------------------
        % 4-PT STENCIL STORAGE ALLOCATION
        % (stencil values in x and y will be stored in these arrays)
        %-----------------------------------------------------------
        deltaX_LagPt_n=zeros(1,4);
        deltaY_LagPt_n=zeros(1,4);

        %-----------------------------------------------------------
        % LOOP ACROSS EACH OF 4 STENCIL POINTS
        %     --> this loop is only necessary for 4-pt Delta 
        %         function from (Peskin 2002)
        %-----------------------------------------------------------
        for nn=1:4

            %-----------------------------------------------------------
            % 'Radii' distances btwn Lag Point and Eulerian Grid Cells
            %-----------------------------------------------------------
            rX=abs(indsX(nn)-xAux);
            rY=abs(indsY(nn)-yAux);

            %-----------------------------------------------------------
            % Evaluate Chosen Delta Function at 'radius' rX or rY
            %    --> Note: (1) divides by dx (or dy) 
            %              (2) default is Peskin 4-pt Delta function 
            %                  (as described in Peskin 2002)
            %-----------------------------------------------------------
            deltaX_LagPt_n(nn) = choose_and_Evaluate_Delta_Function(rX,dx);
            deltaY_LagPt_n(nn) = choose_and_Evaluate_Delta_Function(rY,dy);

        end


        %-------------------------------------------------------
        % Adjust indices to account for periodic boundaries
        %-------------------------------------------------------
        indsX = mod(indsX-1,Nx) +1;
        indsY = mod(indsY-1,Ny) +1;


        %--------------------------------------------------------
        % Get Grid Meshes:
        %      ex) a = [a1 a2 a3]; b = [b1 b2 b3];
        %              [X,Y] = meshgrid(a,b)
        %             _        _          _        _
        %            | a1 a2 a3 |        | b1 b1 b1 |
        %        X = | a1 a2 a3 |    Y = | b2 b2 b2 |        
        %            |_a1 a2 a3_|        |_b3 b3 b3_|
        %
        %   NOTE: Eulerian data structures are handled in IB2d 
        %         using u(j,i) w/ i-x coordinate, j-y-coordinate 
        %--------------------------------------------------------
        %
        % All Eulerian indices in which grid is evaluated for for Lag. Point n
        [indsX,indsY]=meshgrid(indsX,indsY);                    
        %
        % Delta function values for those grid points.
        [deltaX_LagPt_n,deltaY_LagPt_n]=meshgrid(deltaX_LagPt_n,deltaY_LagPt_n); 


        %----------------------------------------------------------------------
        % sub2ind --> takes a matrix index like (5,2) and 
        %             converts it to a single number index
        %                 ex) If Matrix M is R^{10x3}
        %                     M(3,2) = 0.0005  <<--->> M(13) = 0.0005
        %   
        %         --> EulerianIDs gives vector of SINGLE indices to call  
        %             entries out of a matrix that correspond to points in 
        %             Eulerian grid that are used for delta function stencil
        %
        %         --> helps bookkeeping for indices for 'flattening' 2D-array 
        %             into a 1D-column vector 
        %----------------------------------------------------------------------
        EulerianIDs = sub2ind([Ny,Nx],indsY(:),indsX(:));


        %----------------------------------------------------------------
        % Each column in S corresponds to a different lag pt.
        %       --> Each column is of size Nx*Ny and represents
        %           delta function values over ENTIRE (Nx,Ny) grid!
        %       --> deltaX_LagPt_n(:).*deltaY_LagPt_n(:) 
        %           ( product gives the Delta_x * Delta_y values
        %           together over non-zero region of the delta stencil)
        %----------------------------------------------------------------
        S(EulerianIDs,n) = deltaX_LagPt_n(:).*deltaY_LagPt_n(:);


    end  % Ends loop over all Lagrangian Points
  
  
  
  
  
