close all;
clear all;
L = 8;                              % length of computational domain (m)
N = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
dx = L/N;                           % Cartesian mesh width (m)
ds = dx/2;
 
 
 
a=.5;
b=.75;
d=-0.25;
factor_a=.8;
 
F=5e0;
 
theta=zeros(1000,1);
theta_lim=asin(d/b);
theta_test=pi/2;
 
x_points=zeros(1000,1);
z_points=zeros(1000,1);
id_points=zeros(1000,1);
offset = 0;
 
kappa_spring = 1.0e5;               % spring constant (Newton)
kappa_beam = 5.0e3;                % beam stiffness constant (Newton m^2)
%kappa_beam_flexible = kappa_beam/5;               % beam stiffness constant (Newton m^2)
kappa_target = kappa_spring;        % target point penalty spring constant (Newton)
 
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
 
 
 
mesh_name = 'plate2d_whole_';      % structure name
 
xShift = 1;
yShift = 2;
 
x_points=x_points(1:npts)+xShift;
z_points=z_points(1:npts)+yShift;
it_points=id_points(1:npts);
 
plot(x_points(:),z_points(:),'*')
axis([0 8 0 8])

vertex_fid = fopen([mesh_name num2str(N) '.vertex'], 'w');
 
fprintf(vertex_fid, '%d\n', npts);
 
for j=1:npts
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_points(j), z_points(j));
end
 
fclose(vertex_fid);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Write out the spring information for the bell
spring_fid = fopen([mesh_name num2str(N) '.spring'], 'w');
npts_spring_type1=npts-1;
 
fprintf(spring_fid, '%d\n', npts-1);
 
 
for s = 1:c-1
    resting=sqrt((x_points(s)-x_points(s+1))^(2)+(z_points(s)-z_points(s+1))^(2));
    fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_points(s)+1, id_points(s+1)+1, kappa_spring*ds/(ds^2), resting, 0);
end
for s = c+1:npts-1
    resting=sqrt((x_points(s)-x_points(s+1))^(2)+(z_points(s)-z_points(s+1))^(2));
    fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_points(s)+1, id_points(s+1)+1, kappa_spring*ds/(ds^2), resting, 0);
end
resting=sqrt((x_points(1)-x_points(c+1))^(2)+(z_points(1)-z_points(c+1))^(2));
fprintf(spring_fid, '%d %d %1.16e %1.16e %d\n', id_points(1)+1, id_points(c+1)+1, kappa_spring*ds/(ds^2), resting, 0);
 
 
fclose(spring_fid);
 
 
 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Write out the beam information
beam_fid = fopen([mesh_name num2str(N) '.beam'], 'w');
 
fprintf(beam_fid, '%d\n', npts-2);
 
for s = 2:c-1
    C1 = x_points(s-1)+x_points(s+1)-2*x_points(s);
    C2 = z_points(s-1)+z_points(s+1)-2*z_points(s);
    fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(s-1)+1, id_points(s)+1, id_points(s+1)+1, kappa_beam*ds/(ds^4), C1, C2);
end
for s = c+2:npts-1
    C1 = x_points(s-1)+x_points(s+1)-2*x_points(s);
    C2 = z_points(s-1)+z_points(s+1)-2*z_points(s);
    fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(s-1)+1, id_points(s)+1, id_points(s+1)+1, kappa_beam*ds/(ds^4), C1, C2);
end
 
C1 = x_points(c+2)+x_points(1)-2*x_points(c+1);
C2 = z_points(c+2)+z_points(1)-2*z_points(c+1);
fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(c+2)+1, id_points(c+1)+1, id_points(1)+1, kappa_beam*ds/(ds^4), C1, C2);
 
C1 = x_points(c+1)+x_points(2)-2*x_points(1);
C2 = z_points(c+1)+z_points(2)-2*z_points(1);
fprintf(beam_fid, '%d %d %d %1.16e %1.16e %1.16e\n', id_points(c+1)+1, id_points(1)+1, id_points(2)+1, kappa_beam*ds/(ds^4), C1, C2);
 
 
 
fclose(beam_fid);
 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% Write out the muscle information
 
muscle_name= 'plate2d_muscle_';
 
musc_vertex_fid = fopen([muscle_name num2str(N) '.vertex'], 'w');
 
 
hold on
fprintf(musc_vertex_fid, '%d\n', npts_musc*2);
for s = 1:npts_musc
    fprintf(musc_vertex_fid, '%1.16e %1.16e\n', x_points(npts_wing+1-npts_musc+s), z_points(npts_wing+1-npts_musc+s));
    plot(x_points(npts_wing+1-npts_musc+s),z_points(npts_wing+1-npts_musc+s),'r*');
end
for s = 1:npts_musc
    fprintf(musc_vertex_fid, '%1.16e %1.16e\n', x_points(npts-npts_musc+s), z_points(npts-npts_musc+s));
    plot(x_points(npts-npts_musc+s),z_points(npts-npts_musc+s),'r*');
end
fclose(musc_vertex_fid);
 
musc_spring_fid = fopen([muscle_name num2str(N) '.spring'], 'w');
 
fprintf(musc_spring_fid, '%d\n', npts_musc);
 
 
for s = 1:npts_musc
 
    fprintf(musc_spring_fid, '%d %d %1.16e %1.16e %d\n', s-1+1, s+npts_musc-1+1, F, 0, 1);
     
end
 
fclose(musc_spring_fid);