function test_Letters()

xCenter=0.5;
yCenter=0;
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
[xt,yt] = give_Me_The_Letter_Please(ds,len,'t',xCenter-1.05,yCenter);
[xh,yh] = give_Me_The_Letter_Please(ds,len,'h',xCenter-1.2,yCenter);
[xn,yn] = give_Me_The_Letter_Please(ds,len,'n',xCenter-1.375,yCenter);
[xr,yr] = give_Me_The_Letter_Please(ds,len,'r',xCenter-1.575,yCenter);
[xg,yg] = give_Me_The_Letter_Please(ds,len,'g',xCenter-1.775,yCenter);
[xe,ye] = give_Me_The_Letter_Please(ds,len,'e',xCenter-1.985,yCenter);
[xm,ym] = give_Me_The_Letter_Please(ds,len,'m',xCenter-2.25,yCenter);
[xK,yK] = give_Me_The_Letter_Please(ds,len,'K',xCenter-2.55,yCenter);
[xk,yk] = give_Me_The_Letter_Please(ds,len,'k',xCenter-2.75,yCenter);
[xEX,yEX] = give_Me_The_Letter_Please(ds,len,'!',xCenter-2.9,yCenter);
[xQU,yQU] = give_Me_The_Letter_Please(ds,len,'?',xCenter-3.05,yCenter);
[xv,yv] = give_Me_The_Letter_Please(ds,len,'v',xCenter-3.25,yCenter);
[xw,yw] = give_Me_The_Letter_Please(ds,len,'w',xCenter-3.5,yCenter);
[xW,yW] = give_Me_The_Letter_Please(ds,len,'W',xCenter-3.85,yCenter);
[xy,yy] = give_Me_The_Letter_Please(ds,len,'y',xCenter-4.15,yCenter);



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
plot(xK,yK,'r*'); hold on;
plot(xk,yk,'b*'); hold on;
plot(xEX,yEX,'k*'); hold on;
plot(xQU,yQU,'r*'); hold on;
plot(xv,yv,'b*'); hold on;
plot(xw,yw,'k*'); hold on;
plot(xW,yW,'r*'); hold on;
plot(xy,yy,'b*'); hold on;



axis([-1 4 -1 4]);

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
    

elseif strcmp(letter,'k')
    
   yLine = -len/2:ds:len/2;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) -xC-len/3.9;

   XD = 0:ds:len/sqrt(2);
   YD = zeros(1,length(XD));
   
   %rotate!
   for i=1:length(XD)
      xDt(i) = XD(i)*cos(pi/4)-YD(i)*sin(pi/4);
      yDt(i) = XD(i)*sin(pi/4)+YD(i)*cos(pi/4);
      
      xDb(i) = XD(i)*cos(pi/4)-YD(i)*sin(-pi/4);
      yDb(i) = XD(i)*sin(-pi/4)+YD(i)*cos(pi/4);
   end
   xDt = xDt - xC - len/4;
   yDt = yDt - yC;
   
   xDb = xDb - xC - len/4;
   yDb = yDb - yC;
   
   x = [xLine(1:ceil(0.75*end)) xDt(1:ceil(2*end/3)) xDb(1:ceil(2*end/3))];
   y = [yLine(1:ceil(0.75*end)) yDt(1:ceil(2*end/3))-len/8 yDb(1:ceil(2*end/3))-len/5.5];
    
elseif strcmp(letter,'K')
    
   yLine = -len/2:ds:len/2;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) -xC-len/3.9;

   XD = 0:ds:len/sqrt(2);
   YD = zeros(1,length(XD));
   
   %rotate!
   for i=1:length(XD)
      xDt(i) = XD(i)*cos(pi/4)-YD(i)*sin(pi/4);
      yDt(i) = XD(i)*sin(pi/4)+YD(i)*cos(pi/4);
      
      xDb(i) = XD(i)*cos(pi/4)-YD(i)*sin(-pi/4);
      yDb(i) = XD(i)*sin(-pi/4)+YD(i)*cos(pi/4);
   end
   xDt = xDt - xC - len/4;
   yDt = yDt - yC;
   
   xDb = xDb - xC - len/4;
   yDb = yDb - yC;
   
   x = [xLine xDt xDb];
   y = [yLine yDt yDb];

 
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
   xLine = zeros(1,length(yLine));

   xH = -len/4:ds:len/4;
   xH = xH-xC;
   yH = zeros(1,length(xH)) - yC + len/9;
   
   x = [xLine-xC xH];
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

elseif strcmp(letter,'v')
    
   XD = 0:ds:sqrt(5)*len/4;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
   end
   xDr = xDr - xC;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC;
   yDL = yDL - yC - len/2;
 
   x = [xDL xDr(2:end)];
   y = [yDL yDr(2:end)];
   
elseif strcmp(letter,'w')
    
   XD = 0:ds:sqrt(5)*len/4;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
   end
   xDr = xDr - xC - len/4;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC - len/4;
   yDL = yDL - yC - len/2;
 
   x = [xDL xDr(2:end-1) xDL+len/2 xDr(2:end)+len/2];
   y = [yDL yDr(2:end-1) yDL yDr(2:end)];   
   
elseif strcmp(letter,'W')
    
   XD = 0:ds:sqrt(5)*len/2;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
   end
   xDr = xDr - xC - len/4;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC - len/4;
   yDL = yDL - yC - len/2;
 
   x = [xDL xDr(2:floor(end/2)-1) xDL(1:floor(end/2))+len/2 xDr(2:end)+len/2];
   y = [yDL yDr(2:floor(end/2)-1) yDL(1:floor(end/2))       yDr(2:end)];     
   
elseif strcmp(letter,'y')
    
   XD = 0:ds:sqrt(5)*len/4;
   YD = zeros(1,length(XD));
   
   %rotate!
   ang = atan(2);
   for i=1:length(XD)
      xDr(i) = XD(i)*cos(ang)-YD(i)*sin(ang);
      yDr(i) = XD(i)*sin(ang)+YD(i)*cos(ang);
      
      xDL(i) = XD(i)*cos(pi-ang)-YD(i)*sin(pi-ang);
      yDL(i) = XD(i)*sin(pi-ang)+YD(i)*cos(pi-ang);
      
      xB(i) = XD(i)*cos(pi+ang)-YD(i)*sin(pi+ang);
      yB(i) = XD(i)*sin(pi+ang)+YD(i)*cos(pi+ang);
   end
   xDr = xDr - xC;
   yDr = yDr - yC - len/2;
 
   xDL = xDL - xC;
   yDL = yDL - yC - len/2;
   
   xB = xB - xC;
   yB = yB - yC - len/2;
 
   x = [xDL xDr(2:end) xB(2:end)];
   y = [yDL yDr(2:end) yB(2:end)];
      
   
elseif strcmp(letter,'!')
    
   yLine = -len/4:ds:len/2;
   yLine = yLine - yC;
   xLine = zeros(1,length(yLine)) - xC; 
   
    r = len/21;
    dt = ds/(2*r);
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xDot(i) = r*cos(theta(i)) - xC;
       yDot(i) = r*sin(theta(i)) - yC - len/2.1;
    end
    
    x = [xLine xDot];
    y = [yLine yDot];  
    
elseif strcmp(letter,'?')
    
    %xT = 0:ds:len/8-ds; xT = xT - xC - len/4;
    %yT = zeros(1,length(xT)) - yC + len/2;

    r = len/4;
    dt = ds/r;
    theta = 1.25*pi/2:-dt:-pi/2;
    for i=1:length(theta)
       xL(i) = 1.5*r*cos(theta(i)) - xC - len/8;
       yL(i) = r*sin(theta(i)) - yC + len/4;
    end
    
    yV = -ds:-ds:-len/3.25;
    xV = zeros(1,length(yV)) - xC - len/8;
    
    r = len/21;
    dt = ds/(2*r);
    theta = 0:dt:2*pi;
    for i=1:length(theta)
       xDot(i) = r*cos(theta(i)) - xC -len/8;
       yDot(i) = r*sin(theta(i)) - yC - len/2.1;
    end
    
    x = [xL xV xDot];
    y = [yL yV yDot];  
 
end