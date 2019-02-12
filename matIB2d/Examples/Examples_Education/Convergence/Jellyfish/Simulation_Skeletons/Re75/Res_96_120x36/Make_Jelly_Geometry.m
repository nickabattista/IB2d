%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jellyfish Example Courtesy of Alexander P. Hoover, PhD
%
% Converted from IBAMR: 1/16/2018 by NAB.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Make_Jelly_Geometry()

close all;
clear all;
L = 8;                              % height of computational domain (m) for keeping desired resolution
Lh = 10;                            % actual height of computational domain (m) (MATCHES INPUT2D)
Lw = 3;                             % width of computational domain (m) (MATCHES INPUT2D)
N = 96;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = dx/2;
 
a=.5;                               % bell radius (semi-minor axis, horizontal axis, note width=2a)
b=.75;                              % bell semi-major axis 
d=-0.25;
factor_a=.8;
 
F=1e5; %5e0
 
theta=zeros(1000,1);
theta_lim=asin(d/b);
theta_test=pi/2;
 
x_points=zeros(1000,1);
z_points=zeros(1000,1);
id_points=zeros(1000,1);
offset = 0;
 
kappa_spring = 1e7; %1e5               % spring constant (Newton)
kappa_beam = 2.5e5; %1e5    %5e3              % beam stiffness constant (Newton m^2)
%kappa_beam_flexible = kappa_beam/5;   % beam stiffness constant (Newton m^2)
kappa_target = kappa_spring;           % target point penalty spring constant (Newton)
 
c=0;
while(theta_test<(pi-theta_lim))
    c=c+1;
    theta(c)=theta_test;
     
    x_points(c)=a*cos(theta(c));
    z_points(c)=b*sin(theta(c));
    id_points(c)=c-1;
     
    theta_test=ds/((a*sin(theta(c)))^(2)+(b*cos(theta(c)))^(2))^(.5)+theta(c);
     
end
 
c_stiff=c;
 
 
npts=2*c-1;
npts_wing=floor(npts/2);
npts_musc=floor(npts_wing/4);
 
for j=(c+1):(npts)
    x_points(j)=-1*x_points(j-c+1);
    z_points(j)=z_points(j-c+1);
    id_points(j)=j-1;
end
 

mesh_name = 'jelly';
xShift = 1.5;
yShift = 2;
 
x_points=x_points(1:npts)+xShift;
z_points=z_points(1:npts)+yShift;
it_points=id_points(1:npts);
 
plot(x_points(:),z_points(:),'*'); hold on;
axis([0 8 0 8])

% Lag Pts to Mess up Flow At Edge
xBlock = ds:4*ds:Lw-ds;
yBlock = (Lh-5*ds)*ones(1,length(xBlock))+ds;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print .vertex information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');
 
    fprintf(vertex_fid, '%d\n', npts + npts_musc*2 + length(xBlock));
    lag_ct = 0;
    
    %
    % bell
    %
    for j=1:npts
        fprintf(vertex_fid, '%1.16e %1.16e\n', x_points(j), z_points(j));
        lag_ct = lag_ct + 1;
    end
 
    %
    % muscles
    %
    for s = 1:npts_musc
        fprintf(vertex_fid, '%1.16e %1.16e\n', x_points(npts_wing+1-npts_musc+s), z_points(npts_wing+1-npts_musc+s));
        plot(x_points(npts_wing+1-npts_musc+s),z_points(npts_wing+1-npts_musc+s),'r*'); hold on;
        lag_ct = lag_ct + 1;
    end
    for s = 1:npts_musc
        fprintf(vertex_fid, '%1.16e %1.16e\n', x_points(npts-npts_musc+s), z_points(npts-npts_musc+s));
        plot(x_points(npts-npts_musc+s),z_points(npts-npts_musc+s),'r*'); hold on;
        lag_ct = lag_ct + 1;
    end
    
    for ii=1:length(xBlock)
        fprintf(vertex_fid, '%1.16e %1.16e\n', xBlock(ii), yBlock(ii));
    end

fclose(vertex_fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print .spring information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spring_fid = fopen([mesh_name num2str(N) '.spring'], 'w');
    
    npts_spring_type1=npts-1;
 
    fprintf(spring_fid, '%d\n', npts-1 + npts_musc);
 
    fprintf('\nNumber of springs before muscles: %d \n\n',npts-1)
    
    factor = 1;%ds^2/ds;
    
    %
    % bell
    %
    for s = 1:c-1
        resting=sqrt((x_points(s)-x_points(s+1))^(2)+(z_points(s)-z_points(s+1))^(2));
        fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_points(s)+1, id_points(s+1)+1, kappa_spring*ds/(ds^2)*factor, resting, 1);
    end
    for s = c+1:npts-1
        resting=sqrt((x_points(s)-x_points(s+1))^(2)+(z_points(s)-z_points(s+1))^(2));
        fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_points(s)+1, id_points(s+1)+1, kappa_spring*ds/(ds^2)*factor, resting, 1);
    end
    resting=sqrt((x_points(1)-x_points(c+1))^(2)+(z_points(1)-z_points(c+1))^(2));
    fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_points(1)+1, id_points(c+1)+1, kappa_spring*ds/(ds^2)*factor, resting, 1);

    %
    % muscles
    %
    for s = 1:npts_musc
        fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n',npts+s-1+1, npts+s+npts_musc-1+1, F, 0, 1);
    end
 
 
    fclose(spring_fid);
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Print .nonInv_beam information
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beam_fid = fopen([mesh_name num2str(N) '.nonInv_beam'], 'w');
 
    fprintf(beam_fid, '%d\n', npts-2);

    factor=1;% = (ds^4)/ds;
    
    for s = 2:c-1
        C1 = x_points(s-1)+x_points(s+1)-2*x_points(s);
        C2 = z_points(s-1)+z_points(s+1)-2*z_points(s);
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(s-1)+1, id_points(s)+1, id_points(s+1)+1, kappa_beam*ds/(ds^4)*factor, C1, C2);
    end
    for s = c+2:npts-1
        C1 = x_points(s-1)+x_points(s+1)-2*x_points(s);
        C2 = z_points(s-1)+z_points(s+1)-2*z_points(s);
        fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(s-1)+1, id_points(s)+1, id_points(s+1)+1, kappa_beam*ds/(ds^4)*factor, C1, C2);
    end

    C1 = x_points(c+2)+x_points(1)-2*x_points(c+1);
    C2 = z_points(c+2)+z_points(1)-2*z_points(c+1);
    fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(c+2)+1, id_points(c+1)+1, id_points(1)+1, kappa_beam*ds/(ds^4)*factor, C1, C2);

    C1 = x_points(c+1)+x_points(2)-2*x_points(1);
    C2 = z_points(c+1)+z_points(2)-2*z_points(1);
    fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(c+1)+1, id_points(1)+1, id_points(2)+1, kappa_beam*ds/(ds^4)*factor, C1, C2);

 
    fclose(beam_fid);
 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%
% PRINT TARGET POINTS!!!
%
% print target points (flow blocker along edge)
k_Target = 2.5e6;
nBefore = lag_ct; % Counts pts in jellyfish for bookkeeping for .target file
struct_name = ['jelly' num2str(N)];
print_Lagrangian_Target_Pts(xBlock,k_Target,struct_name,nBefore)    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints TARGET points to a file called 'struct_name'.target
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Target_Pts(xLag,k_Target,struct_name,nBefore)

    N = length(xLag);
    Nstart = nBefore+1;
    Nend = nBefore+N;

    target_fid = fopen([struct_name '.target'], 'w');

    fprintf(target_fid, '%d\n', N );

    %Loops over all Lagrangian Pts.
    for s = Nstart:Nend
        fprintf(target_fid, '%d %1.16e\n', s, k_Target);
    end

    fclose(target_fid); 