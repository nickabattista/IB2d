
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: finds coefficients of cubic spline interpolant
%
% INPUTs: (p1,p2): interpolating mediary points
%
% Author: Nick Battista
% Date: February 23, 2018
% Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function interp_Function_Coeffs(p1,p2)

%
% Set up MATRIX for linear system of cubic interpolant
%
mat1 = [1 0 0 0 0 0 0 0 0 0 0 0];
mat2 = [0 1 0 0 0 0 0 0 0 0 0 0];
mat3 = [0 0 1 0 0 0 0 0 0 0 0 0];
mat4 = [1 p1 p1^2 p1^3 -1 -p1 -p1^2 -p1^3 0 0 0 0];

mat5 = [0 1 2*p1 3*p1^2 0 -1 -2*p1 -3*p1^2 0 0 0 0];
mat6 = [0 0 2 6*p1 0 0 -2 -6*p1 0 0 0 0];
mat7 = [0 0 0 0 1 p2 p2^2 p2^3 -1 -p2 -p2^2 -p2^3];
mat8 = [0 0 0 0 0 1 2*p2 3*p2^2 0 -1 -2*p2 -3*p2^2];

mat9 = [0 0 0 0 0 0 2 6*p2 0 0 -2 -6*p2];
mat10= [0 0 0 0 0 0 0 0 1 1 1 1];
mat11= [0 0 0 0 0 0 0 0 0 1 2 3];
mat12= [0 0 0 0 0 0 0 0 0 0 2 6];

mat = [mat1; mat2; mat3; mat4; mat5; mat6; mat7; mat8; mat9; mat10; mat11; mat12];

%
% RHS of Linear System for cubic interpolant
%
rhs = [0 0 0 0 0 0 0 0 0 1 0 0]';

%
% Give coefficients
%
coeffs = mat\rhs


%------------------------------------
% Plot Interpolation Polynomial
%               (and its derivatives)
%------------------------------------

%Plotting t-Values
ds = 0.005;
t = 0:ds:1;
t2=t/2;

%Plot Attributes
ms = 20;

%
% PLOT: interpolating function, g(t)
figure(1)
for i=1:length(t)
    plot(t2(i),g(coeffs,t(i),p1,p2),'.','MarkerSize',ms); hold on;
end
xlabel('t');
ylabel('g(t)');
set(gca,'FontSize',18)

%
% PLOT: interpolating function's derivative, dg(t)/dt
figure(2)
for i=1:length(t)
    plot(t2(i),gP(coeffs,t(i),p1,p2),'.','MarkerSize',ms); hold on;
end
xlabel('t');
ylabel('dg/dt');
set(gca,'FontSize',18)

%
% PLOT: interpolating function's 2nd derivative, d^2g(t)/dt^2
figure(3)
for i=1:length(t)
    plot(t2(i),gPP(coeffs,t(i),p1,p2),'.','MarkerSize',ms); hold on;
end
xlabel('t');
ylabel('d^2g/dt^2');
set(gca,'FontSize',18)



%
% Store coefficients for plotting
%
% vec = []; 
% for i=1:12
%     vec = [vec coeffs(i)];
% end

%strName = 'coeffs_Vec';
%print_Matrix_To_Txt_File(coeffs,strName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION prints matrix to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Matrix_To_Txt_File(a,strName)

nameTxt = [strName '.txt'];

fid = fopen(nameTxt, 'wt'); % Open for writing
for i=1:size(a,1)
   fprintf(fid, '%.14f ', a(i,:));
   fprintf(fid, '\n');
end
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluate interpolating polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = g(coeffs,x,p1,p2)

if x<=p1
    val = coeffs(1) + coeffs(2)*x + coeffs(3)*x^2 + coeffs(4)*x^3;
elseif x<=p2
    val = coeffs(5) + coeffs(6)*x + coeffs(7)*x^2 + coeffs(8)*x^3;
else
    val = coeffs(9) + coeffs(10)*x + coeffs(11)*x^2 + coeffs(12)*x^3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluate interpolating polynomial 1st derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = gP(coeffs,x,p1,p2)

if x<=p1
    val = coeffs(2) + 2*coeffs(3)*x + 3*coeffs(4)*x^2;
elseif x<=p2
    val = coeffs(6) + 2*coeffs(7)*x + 3*coeffs(8)*x^2;
else
    val = coeffs(10) + 2*coeffs(11)*x + 3*coeffs(12)*x^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluate interpolating polynomial 2nd derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = gPP(coeffs,x,p1,p2)

if x<=p1
    val = 2*coeffs(3) + 6*coeffs(4)*x;
elseif x<=p2
    val = 2*coeffs(7) + 6*coeffs(8)*x;
else
    val = 2*coeffs(11) + 6*coeffs(12)*x;
end

