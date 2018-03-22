function interp_Function_Coeffs(x1,x2)

mat1 = [-x1^2 0 x1^3 x1^2 x1 1];
mat2 = [-2*x1 0 3*x1^2 2*x1 1 0];
mat3 = [-2 0 6*x1 2 0 0];
mat4 = [0 (x2-1)^2 x2^3 x2^2 x2 1];
mat5 = [0 2*(x2-1) 3*x2^2 2*x2 1 0];
mat6 = [0 2 6*x2 2 0 0];

mat = [mat1; mat2; mat3; mat4; mat5; mat6];

rhs = [0 0 0 1 0 0]';

coeffs = mat\rhs

a = coeffs(1);
b = coeffs(2); 
c = coeffs(3);
d = coeffs(4);
g = coeffs(5); 
h = coeffs(6);

ds = 0.01;
x1_part = 0:ds:x1;
x2_part = x1:ds:x2;
x3_part = x2:ds:1.0;

%
% PLOT INTERPOLATING FUNCTION
%
figure(1)
for i=1:length(x1_part)
    x = x1_part(i);
    plot(x,a*x^2,'ro'); hold on;
end

for i=1:length(x2_part)
    x = x2_part(i);
    plot(x,c*x^3+d*x^2+g*x+h,'*'); hold on;
end

for i=1:length(x3_part)
   x = x3_part(i);
   plot(x,-b*(x-1)^2+1,'go'); hold on;
end

xtotal = 0:ds:1;
for i=1:length(xtotal)
   x = xtotal(i);
   plot(x,0.5*(1+tanh(4.3*(x-0.5))),'m*'); hold on;
end
title('function vals');


%
% PLOT 1ST DERIVATIVE OF INTERPOLATING FUNCTION
%
figure(2)
for i=1:length(x1_part)
    x = x1_part(i);
    plot(x,2*a*x,'ro'); hold on;
end

for i=1:length(x2_part)
    x = x2_part(i);
    plot(x,3*c*x^2+2*d*x+g,'*'); hold on;
end

for i=1:length(x3_part)
   x = x3_part(i);
   plot(x,-2*b*(x-1),'go'); hold on;
end
title('1st deriv');


%
% PLOT2ND DERIVATIVE OF INTERPOLATING FUNCTION
%
figure(3) 
for i=1:length(x1_part)
    x = x1_part(i);
    plot(x,2*a,'ro'); hold on;
end

for i=1:length(x2_part)
    x = x2_part(i);
    plot(x,6*c*x+2*d,'*'); hold on;
end

for i=1:length(x3_part)
   x = x3_part(i);
   plot(x,-2*b,'go'); hold on;
end
title('2nd deriv');
