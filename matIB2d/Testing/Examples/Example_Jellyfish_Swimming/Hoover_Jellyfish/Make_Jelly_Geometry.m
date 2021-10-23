%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jellyfish Example Courtesy of Alexander P. Hoover, PhD and 
%                                           Laura A. Miller, PhD
%
%       -> Converted from IBAMR: 1/16/2018 by NAB.
%
%       -> Modified from original on 3/5/2021 by NAB.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Make_Jelly_Geometry()


close all;

%-----------------------------------------------
% Grid Parameters
%-----------------------------------------------
Ly = 10;                  % height of computational domain (m) (MATCHES INPUT2D)
Lx = 3;                  % width of computational domain (m) (MATCHES INPUT2D)
Ny = 640;                % number of Cartesian grid meshwidths 
dx = Ly/Ny;              % Cartesian mesh resolution (width(m))
ds = dx/2;               % Lagrangian Pt Spacing



%-----------------------------------------------
% Jellyfish Geometric Parameters (ellipse)
%-----------------------------------------------
a=.5;                     % bell radius (semi-minor axis, horizontal axis, note width=2a)
b=.75;                    % bell semi-major axis 
d=-0.25;                  % lowest point along bell
 
 

%-----------------------------------------------
% Set up theta values for creating jellyfish
%-----------------------------------------------
theta_lim=asin(d/b);  % last theta value for jellyfish (on left side)
theta_test=pi/2;      % starting theta value (top of jellyfish)
 

%-----------------------------------------------
% Spring and Beam Stiffness Coefficients; 
%                       Jellyfish Muscle Force
%-----------------------------------------------
kappa_spring = 1e7;       % spring constant (Newton)
kappa_beam = 2.5e5;       % beam stiffness constant (Newton m^2)
F=1e5;                    % JELLYFISH CONTRACTION FORCE


%-----------------------------------------------
% CREATE **LEFT** SIDE JELLYFISH GEOMETRY!
%-----------------------------------------------
c=0;
while(theta_test<(pi-theta_lim))
    c=c+1;
    theta(c)=theta_test;
     
    xL(c)=a*cos(theta(c));
    yL(c)=b*sin(theta(c));
    id_points(c)=c-1;
     
    theta_test=ds/((a*sin(theta(c)))^(2)+(b*cos(theta(c)))^(2))^(.5)+theta(c);
     
end


%-----------------------------------------------
% NUMBER OF PTS ALONG JELLYFISH AND MUSCLES
%-----------------------------------------------
npts= 2*length(xL)-1;        % # OF POINTS ALONG ENTIRE JELLYFISH
npts_wing=length(xL)-1;      % # OF PTS ALONG ONE SIDE OF JELLYFISH (NOT INCLUDING CENTER PT)   
npts_musc=floor(npts_wing/4);% # OF MUSCLES TO CONTRACT BELL
 

%-----------------------------------------------
% CREATE RIGHT SIDE JELLYFISH GEOMETRY!
%-----------------------------------------------
xR = -xL(2:end);  % So not to include the CENTER point
yR = yL(2:end);   % So not to include the CENTER point


%-----------------------------------------------
% TRANSLATE JELLYFISH
%-----------------------------------------------
xShift = 1.5;
yShift = 2.0; 
xLag=[xL xR]+xShift;
yLag=[yL yR]+yShift;

  

%-----------------------------------------------
% STARTING INDEX FOR MUSCLES ON LEFT AND RIGHT
%-----------------------------------------------
ind_Musc_L = 1 + npts_wing - npts_musc + 1;
ind_Musc_R = 1 + 2*npts_wing - npts_musc + 1;





%-----------------------------------------------
% STORE NAME
%-----------------------------------------------
mesh_name = 'jelly';





%-----------------------------------------------
% PLOT/TEST GEOMETRY
%-----------------------------------------------
plot(xLag,yLag,'*'); hold on;
axis([0 8 0 8])


%-----------------------------------------------
% Lag Pts to Occlude Flow At Edge
%-----------------------------------------------
xBlock = ds:4*ds:Lx-ds;
yBlock = (Ly-5*ds)*ones(1,length(xBlock))+ds;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print .vertex information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen([mesh_name '.vertex'], 'w');
 
    fprintf(vertex_fid, '%d\n', npts + length(xBlock) );
    lag_ct = 0;
    
    %
    % JELLYFISH 1'S BELL
    %
    for j=1:npts
        fprintf(vertex_fid, '%1.16e %1.16e\n', xLag(j), yLag(j));
        lag_ct = lag_ct + 1;
    end
    
    
    %
    % flow blocker across top edge
    %
    for ii=1:length(xBlock)
        fprintf(vertex_fid, '%1.16e %1.16e\n', xBlock(ii), yBlock(ii));
    end

fclose(vertex_fid);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print .spring information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spring_fid = fopen([mesh_name '.spring'], 'w');
    
 
    fprintf(spring_fid, '%d\n', npts-1 + npts_musc);
 
    fprintf('\n         ***** For update_Springs.m file ***** \n\n');
    %fprintf('  -> Number of springs BEFORE muscles on Jelly-1: %d \n',npts-1)
    fprintf('  -> Spring-IDs:\n');
    fprintf('                 Jelly-1 Muscles: %d to %d \n\n',npts,npts+npts_musc-1)

        
    factor = 1;%ds^2/ds;
 
    
    %
    % ----------------- JELLYFISH 1 --------------------
    %
    
    %
    % JELLYFISH 1: LEFT SIDE OF BELL
    %
    for s = 1:npts_wing % NOTE: attaches spring along each point on left side
        id_1 = s;
        id_2 = s+1;
        dsRest = sqrt( ( xLag(id_1)-xLag(id_2) )^2 + ( yLag(id_1)-yLag(id_2) )^2 );
        fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_1, id_2, kappa_spring*ds/(ds^2)*factor, dsRest, 1);
    end
    
    %
    % JELLYFISH 1: RIGHT SIDE OF BELL
    %
    for s = 1:npts_wing
        if s==1
            id_1 = 1;
            id_2 = (1+npts_wing)+1; % Note: (1+npts_wing) gives total # of pts on left side)
            dsRest = sqrt( ( xLag(id_1)-xLag(id_2) )^2 + ( yLag(id_1)-yLag(id_2) )^2 );
            fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_1, id_2, kappa_spring*ds/(ds^2)*factor, dsRest, 1);    
        else
            id_1 = (s-1) + (npts_wing+1);
            id_2 = s + (npts_wing+1);
            dsRest = sqrt( ( xLag(id_1)-xLag(id_2) )^2 + ( yLag(id_1)-yLag(id_2) )^2 );
            fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_1, id_2, kappa_spring*ds/(ds^2)*factor, dsRest, 1);
        end
    end
    
    %
    % JELLYFISH 1: MUSCLES between ends of bell
    %
    for s = 0:npts_musc-1
        id_1 = ind_Musc_L+s;
        id_2 = ind_Musc_R+s;
        dsRest = sqrt( ( xLag(id_1)-xLag(id_2) )^2 + ( yLag(id_1)-yLag(id_2) )^2 );
        fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n',id_1, id_2, F, dsRest, 1);
    end


    fclose(spring_fid);
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print .nonInv_beam information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam_fid = fopen([mesh_name '.nonInv_beam'], 'w');
 
    fprintf(beam_fid, '%d\n', npts-2);

    factor=1;% = (ds^4)/ds;
    
    %
    % ----------------- JELLYFISH 1 --------------------
    %
    
    %
    % JELLYFISH 1: BEAM BTWN MIDDLE THREE PTS AT TOP OF BELL
    %
    id_L = 2;                  % First point on left sides
    id_M = 1;                  % Middlemost point of bell
    id_R = (1+npts_wing) + 1;  % First point on right side
    C1 = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
    C2 = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
    fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kappa_beam*ds/(ds^4)*factor, C1, C2);

    %
    % JELLYFISH 1: BEAMS ALONG LEFT SIDE OF BELL
    %
    for s = 2:npts_wing
            id_L = s+1;      % Left Point ID
            id_M = s;        % Middle Point ID
            id_R = s-1;      % Right Point ID
            C1 = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
            C2 = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
            fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kappa_beam*ds/(ds^4)*factor, C1, C2);
    end
    
    
    %
    % JELLYFISH 1: BEAMS ALONG RIGHT SIDE OF BELL
    %
    for s = 2:npts_wing
        
            % auxiliary index for consistency with *LEFT* side
            ss = s + (npts_wing);
        
            if s==2
                id_L = 1;      % MIDDLE POINT
                id_M = ss;     % FIRST POINT ON RIGHT SIDE AFTER MIDDLE PTs
                id_R = ss+1;   % 2ND POINT ON RIGHT SIDE AFTER MIDDLE PT
                C1 = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
                C2 = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kappa_beam*ds/(ds^4)*factor, C1, C2);
            else
                id_L = ss-1;      % Left Point ID
                id_M = ss;        % Middle Point ID
                id_R = ss+1;      % Right Point ID
                C1 = xLag(id_L)+xLag(id_R)-2*xLag(id_M);
                C2 = yLag(id_L)+yLag(id_R)-2*yLag(id_M);  
                fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_L, id_M, id_R, kappa_beam*ds/(ds^4)*factor, C1, C2);    
            end
                    
    end
  
 
    fclose(beam_fid);
 

    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%
% PRINT TARGET POINTS!!!
%
% print target points (flow blocker along edge)
k_Target = 2.5e6;
nBefore = length(xLag); % Counts pts in jellyfish for bookkeeping for .target file
print_Lagrangian_Target_Pts(xBlock,k_Target,mesh_name,nBefore)    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called 'struct_name'.target
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xBlock,k_Target,struct_name,nBefore)

    N = length(xBlock);
    Nstart = nBefore+1;
    Nend = nBefore+N;

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = Nstart:Nend
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid); 
    
