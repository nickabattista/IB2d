function interp_Function_Coeffs(x1,x2)


mat1 = [1 0 0 0 0 0 0 0 0 0 0 0];
mat2 = [0 1 0 0 0 0 0 0 0 0 0 0];
mat3 = [0 0 1 0 0 0 0 0 0 0 0 0];
mat4 = [1 x1 x1^2 x1^3 -1 -x1 -x1^2 -x1^3 0 0 0 0];
mat5 = [0 1 2*x1 3*x1^2 0 -1 -2*x1 -3*x1^2 0 0 0 0];
mat6 = [0 0 2 6*x1 0 0 -2 -6*x1 0 0 0 0];

mat7 = [0 0 0 0 1 x2 x2^2 x2^3 -1 -x2 -x2^2 -x2^3];
mat8 = [0 0 0 0 0 1 2*x2 3*x2^2 0 -1 -2*x2 -3*x2^2];
mat9 = [0 0 0 0 0 0 2 6*x2 0 0 -2 -6*x2];

mat10= [0 0 0 0 0 0 0 0 1 1 1 1];
mat11= [0 0 0 0 0 0 0 0 0 1 2 3];
mat12= [0 0 0 0 0 0 0 0 0 0 2 6];

mat = [mat1; mat2; mat3; mat4; mat5; mat6; mat7; mat8; mat9; mat10; mat11; mat12];

rhs = [0 0 0 0 0 0 0 0 0 1 0 0]';

%mat1 = [-x1^2 0 x1^3 x1^2 x1 1];
%mat2 = [-2*x1 0 3*x1^2 2*x1 1 0];
%mat3 = [-2 0 6*x1 2 0 0];
%mat4 = [0 (x2-1)^2 x2^3 x2^2 x2 1];
%mat5 = [0 2*(x2-1) 3*x2^2 2*x2 1 0];
%mat6 = [0 2 6*x2 2 0 0];
%mat = [mat1; mat2; mat3; mat4; mat5; mat6];
%rhs = [0 0 0 1 0 0]';

coeffs = mat\rhs



a0 = coeffs(1);
a1 = coeffs(2); 
a2 = coeffs(3);
a3 = coeffs(4);
b0 = coeffs(5); 
b1 = coeffs(6);
b2 = coeffs(7);
b3 = coeffs(8); 
c0 = coeffs(9);
c1 = coeffs(10);
c2 = coeffs(11); 
c3 = coeffs(12);

ds = 0.01;
x = 0:ds:1;

figure(1)
for i=1:length(x)
    plot(x(i),g(coeffs,x(i),x1,x2),'*'); hold on;
end


figure(2)
for i=1:length(x)
    plot(x(i),gP(coeffs,x(i),x1,x2),'*'); hold on;
end

figure(3)
for i=1:length(x)
    plot(x(i),gPP(coeffs,x(i),x1,x2),'*'); hold on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluate interpolating polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = g(coeffs,x,x1,x2)

if x<=x1
    val = coeffs(1) + coeffs(2)*x + coeffs(3)*x^2 + coeffs(4)*x^3;
elseif x<=x2
    val = coeffs(5) + coeffs(6)*x + coeffs(7)*x^2 + coeffs(8)*x^3;
else
    val = coeffs(9) + coeffs(10)*x + coeffs(11)*x^2 + coeffs(12)*x^3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluate interpolating polynomial 1st derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = gP(coeffs,x,x1,x2)

if x<=x1
    val = coeffs(2) + 2*coeffs(3)*x + 3*coeffs(4)*x^2;
elseif x<=x2
    val = coeffs(6) + 2*coeffs(7)*x + 3*coeffs(8)*x^2;
else
    val = coeffs(10) + 2*coeffs(11)*x + 3*coeffs(12)*x^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: evaluate interpolating polynomial 2nd derivative
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = gPP(coeffs,x,x1,x2)

if x<=x1
    val = 2*coeffs(3) + 6*coeffs(4)*x;
elseif x<=x2
    val = 2*coeffs(7) + 6*coeffs(8)*x;
else
    val = 2*coeffs(11) + 6*coeffs(12)*x;
end

% 
% figure(1)
% for i=1:length(x1_part)
%     x = x1_part(i);
%     plot(x,a*x^2,'ro'); hold on;
% end
% 
% for i=1:length(x2_part)
%     x = x2_part(i);
%     plot(x,c*x^3+d*x^2+g*x+h,'*'); hold on;
% end
% 
% for i=1:length(x3_part)
%    x = x3_part(i);
%    plot(x,-b*(x-1)^2+1,'go'); hold on;
% end
% 
% xtotal = 0:ds:1;
% for i=1:length(xtotal)
%    x = xtotal(i);
%    plot(x,0.5*(1+tanh(4.3*(x-0.5))),'m*'); hold on;
% end
% title('function vals');
% 
% 
% figure(2)
% for i=1:length(x1_part)
%     x = x1_part(i);
%     plot(x,2*a*x,'ro'); hold on;
% end
% 
% for i=1:length(x2_part)
%     x = x2_part(i);
%     plot(x,3*c*x^2+2*d*x+g,'*'); hold on;
% end
% 
% for i=1:length(x3_part)
%    x = x3_part(i);
%    plot(x,-2*b*(x-1),'go'); hold on;
% end
% title('1st deriv');
% 
% 
% figure(3)
% for i=1:length(x1_part)
%     x = x1_part(i);
%     plot(x,2*a,'ro'); hold on;
% end
% 
% for i=1:length(x2_part)
%     x = x2_part(i);
%     plot(x,6*c*x+2*d,'*'); hold on;
% end
% 
% for i=1:length(x3_part)
%    x = x3_part(i);
%    plot(x,-2*b,'go'); hold on;
% end
% title('2nd deriv');
