function test_Letters()

xCenter=1;
yCenter=0.5;
len=0.25;
ds=0.0125;

[xa,ya] = give_Me_The_Letter_Please(ds,len,'a',xCenter-0.125,yCenter);
[xC,yC] = give_Me_The_Letter_Please(ds,len,'C',xCenter,yCenter);
[xH,yH] = give_Me_The_Letter_Please(ds,len,'H',xCenter-0.3,yCenter);
[xi,yi] = give_Me_The_Letter_Please(ds,len,'i',xCenter-0.425,yCenter);
[xo,yo] = give_Me_The_Letter_Please(ds,len,'o',xCenter-0.525,yCenter);
[xl,yl] = give_Me_The_Letter_Please(ds,len,'l',xCenter-0.7,yCenter);
[xd,yd] = give_Me_The_Letter_Please(ds,len,'d',xCenter-0.75,yCenter);
[xu,yu] = give_Me_The_Letter_Please(ds,len,'u',xCenter-0.925,yCenter);
[xt,yt] = give_Me_The_Letter_Please(ds,len,'t',xCenter-1.15,yCenter);
[xh,yh] = give_Me_The_Letter_Please(ds,len,'h',xCenter-1.3,yCenter);
[xn,yn] = give_Me_The_Letter_Please(ds,len,'n',xCenter-1.475,yCenter);
[xr,yr] = give_Me_The_Letter_Please(ds,len,'r',xCenter-1.675,yCenter);
[xg,yg] = give_Me_The_Letter_Please(ds,len,'g',xCenter-1.875,yCenter);
[xe,ye] = give_Me_The_Letter_Please(ds,len,'e',xCenter-2.085,yCenter);
[xm,ym] = give_Me_The_Letter_Please(ds,len,'m',xCenter-2.35,yCenter);



plot(xa,ya,'r*'); hold on;
plot(xC,yC,'b*'); hold on;
plot(xH,yH,'k*'); hold on;
plot(xi,yi,'r*'); hold on;
plot(xo,yo,'b*'); hold on;
plot(xl,yl,'k*'); hold on;
plot(xd,yd,'r*'); hold on;
plot(xu,yu,'b*'); hold on;
plot(xt,yt,'k*'); hold on;
plot(xh,yh,'r*'); hold on;
plot(xn,yn,'b*'); hold on;
plot(xr,yr,'k*'); hold on;
plot(xg,yg,'r*'); hold on;
plot(xe,ye,'b*'); hold on;
plot(xm,ym,'k*'); hold on;

axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give me the letter!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y] = give_Me_The_Letter_Please(ds,len,letter,xC,yC)

if strcmp(letter,'a')
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xA(i) = r*cos(theta(i)) - xC;
       yA(i) = r*sin(theta(i)) - yC - len/4;
    end
    
    yLine = [-r:ds:r r]; yLine = yLine - yC - len/4;
    xLine = r*ones(1,length(yLine)) - (xC-ds/6);
    
    x = [xA xLine];
    y = [yA yLine];

elseif strcmp(letter,'C')
    
    r = len/2;
        dt = ds/r;
    theta = pi/2:dt:pi-0.25*dt;
    theta = [theta pi -theta];
    for i=1:length(theta)
       x(i) = 3*r/4*cos(theta(i)) - xC;
       y(i) = r*sin(theta(i)) - yC;
    end
    
    xTB = 0:3*ds/4:len/6;
    xTB = xTB - xC;
    yT = r*ones(1,length(xTB))-yC;
    yB = -r*ones(1,length(xTB))-yC;
    x = [xTB x xTB];
    y = [yT y yB];
    
elseif strcmp(letter,'d')
    
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xA(i) = r*cos(theta(i)) - xC;
       yA(i) = r*sin(theta(i)) - yC - len/4;
    end
    
    yLine = [-len/2:ds:len/2.25 r]; yLine = yLine - yC;
    xLine = r*ones(1,length(yLine)) - (xC-ds/6);
    
    x = [xA xLine];
    y = [yA yLine];
    
 elseif strcmp(letter,'e')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC;
    yH2 = yH - len/4;
    yH3 = yH - len/2;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    yV3 = 0:-ds:-len/4;
    yV3 = yV3 - yC;
    xV3 = r*ones(1,length(yV3));
        
    x = [-xV-(xC-ds/6) xH xV3-(xC-ds/6) xH xH];
    y = [yV yH yV3 yH2 yH3];
           
    
elseif strcmp(letter,'g')
    
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xA(i) = r*cos(theta(i)) - xC;
       yA(i) = r*sin(theta(i)) - yC - len/4;
    end
    
    yLine = [-len/2:ds:len/2.25 r]; yLine = yLine - yC - len/3;
    xLine = r*ones(1,length(yLine)) - (xC-ds/6);
    
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC - 0.85*len;
    
    x = [xA xLine xH];
    y = [yA yLine yH];
    
elseif strcmp(letter,'h')

   yLine = -len/2:ds:0-ds/4;
   yLine = [yLine 0 -yLine];
   yLine1 = yLine - yC;
   xLine1 = 1/2*len/2*ones(1,length(yLine1));
   
   yLine2 = -len/2:ds:0-ds/4;
   yLine2 = [yLine2 0];
   yLine2 = yLine2 - yC;
   xLine2 = 1/2*len/2*ones(1,length(yLine2));
   
   xH = (-1/2*len/2+ds/2:ds:1/2*len/2-ds/4);
   xH = xH - xC;
   yH = zeros(1,length(xH)) - yC;
   
   x = [-xLine1-xC xH xLine2-xC];
   y = [yLine1 yH yLine2];    
   
elseif strcmp(letter,'H')
   yLine = -len/2:ds:0-ds/4;
   yLine = [yLine 0 -yLine];
   yLine = yLine - yC;
   xLine = 1/2*len/2*ones(1,length(yLine));
   xH = (-1/2*len/2+ds/2:ds:1/2*len/2-ds/4);
   xH = xH - xC;
   yH = zeros(1,length(xH)) - yC;
   
   x = [-xLine-xC xH xLine-xC];
   y = [yLine yH yLine];
   
elseif strcmp(letter,'i')
    
   yLine = -len/2:ds:0-ds/4;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) - xC; 
   
    r = len/21;
    dt = ds/(2*r);
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xDot(i) = r*cos(theta(i)) - xC;
       yDot(i) = r*sin(theta(i)) - yC + len/7;
    end
    
    x = [xLine xDot];
    y = [yLine yDot];
 
elseif strcmp(letter,'l')
    
   yLine = -len/2:ds:len/2.25;
   yLine = yLine - yC;
   xLine = 1/2*len/2*ones(1,length(yLine));

   x = -xLine-xC;
   y = yLine;

elseif strcmp(letter,'m')
        
    r = len/2;
    xH = -r:ds:r;  xH = xH - xC;
    yH = zeros(1,length(xH)) - yC;
    yV = 0:-ds:-len/2; yV = yV - yC;
    xV = -r*ones(1,length(yV)) - xC;
    xV2 = zeros(1,length(yV)) - xC;
    xV3 = r*ones(1,length(yV)) - xC;
    
    yVS = 0:ds:len/12;
    yVS = yVS - yC;
    xVS = -r*ones(1,length(yVS)) - xC;
        
    x = [xVS xH xV xV2 xV3];
    y = [yVS yH yV yV yV];
    
elseif strcmp(letter,'n')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    yV2 = 0:ds:len/12;
    yV2 = yV2 - yC;
    xV2 = r*ones(1,length(yV2));
        
    x = [-xV-(xC-ds/6) -xV2-(xC-ds/6) xH xV-(xC-ds/6)];
    y = [yV yV2 yH yV];
    
   
elseif strcmp(letter,'o')
       
    r = len/4.1;
    dt = ds/r;
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       x(i) = r*cos(theta(i)) - xC;
       y(i) = r*sin(theta(i)) - yC - len/4;
    end
    
 elseif strcmp(letter,'r')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    yV2 = 0:ds:len/12;
    yV2 = yV2 - yC;
    xV2 = r*ones(1,length(yV2));
    
    yV3 = 0:-ds:-len/10;
    yV3 = yV3 - yC;
    xV3 = r*ones(1,length(yV3));
        
    x = [-xV-(xC-ds/6) -xV2-(xC-ds/6) xH xV3-(xC-ds/6)];
    y = [yV yV2 yH yV3];
       
    
elseif strcmp(letter,'t')
    
   yLine = -len/2:ds:len/2.25;
   yLine = yLine - yC;
   xLine = 1/2*len/2*ones(1,length(yLine));

   xH = -len/6.5:ds:len/6.5 - xC;
   yH = zeros(1,length(xH)) - yC + len/9;
   
   x = [-xLine-xC xH];
   y = [yLine yH];
   
    
elseif strcmp(letter,'u')
        
    r = len/4.1;
    yV = [-r:ds:r r]; 
    xH = yV - xC;
    yH = zeros(1,length(xH)) - yC - len/2;
    
    yV = yV - yC - len/4;
    xV = r*ones(1,length(yV));
    
    x = [-xV-(xC-ds/6) xH xV-(xC-ds/6)];
    y = [yV yH yV];
    
end