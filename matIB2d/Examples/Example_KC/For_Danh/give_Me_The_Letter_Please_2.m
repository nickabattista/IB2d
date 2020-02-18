%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: give me the letter!!!!!
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

elseif strcmp(letter,'J')
    
    xFlat=-len/2:ds:len/2;
    yFlat=len/2*ones(size(xFlat));
    
    L_vert = 2*len/3;
    xShift=len/5;
    yShift=len/2-L_vert;
    yVert=len/2:-ds:yShift;
    %yVert = yVert + yShift;
    xVert=len/5*ones(size(yVert));
    
    r=len/3;
    arcL = pi*r;
    lenVec=-arcL:ds:0;
    N = length(lenVec);
    dAng = pi/N;
    angVec=-pi:dAng:0;
    for i=1:N
               
        xCirc(i) = r*cos( angVec(i) );
        yCirc(i) = r*sin( angVec(i) );
        
    end
    xCirc = xCirc-r+xShift;
    yCirc = yCirc+yShift;
    
    
    x = [xCirc xVert xFlat];
    y = [yCirc yVert yFlat];
    
    x = x+xC;
    y = y+yC;
    
    
elseif strcmp(letter,'U')    

    
    L_vert = 2*len/3;
    xShift=len/2;
    yShift=len/2-L_vert;
    yVert=len/2:-ds:yShift;
    %yVert = yVert + yShift;
    xVert=-len/3*ones(size(yVert));
    xVert2=xVert+2*len/3;
    
    r=len/3;
    arcL = pi*r;
    lenVec=-arcL:ds:0;
    N = length(lenVec);
    dAng = pi/N;
    angVec=-pi:dAng:0;
    for i=1:N
               
        xCirc(i) = r*cos( angVec(i) );
        yCirc(i) = r*sin( angVec(i) );
        
    end
    xCirc = xCirc;%-r%+xShift;
    yCirc = yCirc+yShift;
    
    
    x = [xVert xVert2 xCirc];
    y = [yVert yVert  yCirc];
    
    x = x+xC;
    y = y+yC;
    
    
    
elseif strcmp(letter,'N')
   
    yVert=-len/2:ds:len/2;
    xVert=zeros(size(yVert));
    xVert2=xVert+2/3*len;
    
    x1=xVert(end); y1 = yVert(end);
    x2=xVert2(1);  y2 = yVert(1);
    dist = sqrt( (x1-x2)^2 + (y1-y2)^2 );
    
    xD = 0:ds:dist;
    xD = xD-dist/2;
    yD = zeros(size(xD));
    
    ang=atan(2*y2/x2);
    for i=1:length(xD)       
        xD2(i) = xD(i)*cos( ang ) - yD(i)*sin( ang );
        yD2(i) = xD(i)*sin( ang ) + yD(i)*cos( ang );
    end
    xD2 = xD2 + dist/3.5;
    
    x = [xVert xD2 xVert2];
    y = [yVert yD2 yVert];
    
    x = x+xC;
    y = y+yC;
    
elseif strcmp(letter,'T')
    
    xFlat=-len/2:ds:len/2;
    yFlat=len/2*ones(size(xFlat));
    
    yVert=-len/2:ds:len/2-ds;
    xVert=zeros(size(yVert));
    
    x = [xFlat xVert];
    y = [yFlat yVert];
    
    x = x+xC;
    y = y+yC;
    

elseif strcmp(letter,'C')
    
    r = len/2;
        dt = ds/r;
    theta = pi/2:dt:pi-0.25*dt;
    theta = [theta pi -theta];
    for i=1:length(theta)
       x(i) = 3*r/4*cos(theta(i)) - xC;
       y(i) = r*sin(theta(i)) - yC;
    end
    
    xTB = 0:3*ds/4:len/4;
    xTB = xTB - xC;
    yT = r*ones(1,length(xTB))-yC;
    yB = -r*ones(1,length(xTB))-yC;
    x = [xTB x xTB];
    y = [yT y yB];
    
    x = x+2*xC;
    y = y+2*yC;
    
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
   yLine = yLine; %- yC;
   xLine = zeros(1,length(yLine))-len/3.9; %-xC;

   XD = 0:ds:len/sqrt(2);
   YD = zeros(1,length(XD));
   
   %rotate!
   for i=1:length(XD)
      xDt(i) = XD(i)*cos(pi/4)-YD(i)*sin(pi/4);
      yDt(i) = XD(i)*sin(pi/4)+YD(i)*cos(pi/4);
      
      xDb(i) = XD(i)*cos(pi/4)-YD(i)*sin(-pi/4);
      yDb(i) = XD(i)*sin(-pi/4)+YD(i)*cos(pi/4);
   end
   xDt = xDt  - len/4;
   yDt = yDt ;
   
   xDb = xDb   - len/4;
   yDb = yDb ;
   
   x = [xLine xDt xDb];
   y = [yLine yDt yDb];
   
   x=x+xC;
   y=y+yC;

 
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
 
   x = [xDL xDr(2:end) xDL+len/2 xDr(2:end)+len/2];
   y = [yDL yDr(2:end) yDL yDr(2:end)];   
   
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
 
   x = [xDL xDr(2:floor(end/2)) xDL(1:floor(end/2))+len/2 xDr(2:end)+len/2];
   y = [yDL yDr(2:floor(end/2)) yDL(1:floor(end/2))       yDr(2:end)];     
   
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