%
% SHIFT VALVE GEOMETRY INTO PLACE!
%
for j=1:5
   Ninfo(j,1) = (2*j-1)*length(xE);                     % Last index of point on chamber
   Ninfo(j,2) = lenOld + (j-1)*(length(xV)-2) + 1;   % 1st index of ACTUAL valve
   Ninfo(j,3) = Ninfo(j,2)+( length( xV ) - 2 )/2;      % Starting index of OTHER SIDE Valve
   Ninfo(j,4) = Ninfo(j,2)+length( xV ) - 2;         % Last point along VALVE
   Ninfo(j,5) = 2*j*length(xE);                           % BOTTOM OF CHAMBER
   %
   xLag = [xLag xV(2:lenVhalf)+xE(end)+(j-1)*xOff xV(lenVhalf+2:end)+xE(end)+(j-1)*xOff];
   yLag = [yLag yV(2:lenVhalf)+3.0-d/2 yV(lenVhalf+2:end)+3.0-d/2];
   %
   plot( xLag( Ninfo(j,1) ), yLag( Ninfo(j,1) ), 'r.','MarkerSize',30); hold on;
   plot( xLag( Ninfo(j,5) ), yLag( Ninfo(j,5) ), 'b.','MarkerSize',30); hold on;

end

Ninfo(6,1) = 1;
Ninfo(6,2) = lenOld+1;
Ninfo(6,3) = lenOld +( length( xV ) - 2 )/2;
Ninfo(6,4) = lenOld +length(xV)-2;
Ninfo(6,5) = 1 + length(xE);