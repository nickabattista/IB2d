function plot_Bell_Contraction()

offset = 0.125;
tP = 0.25;
total_time = 3*tP + 2*offset;
dt = tP/125;


%Coefficients for Polynomial Phase-Interpolation
a = 2.739726027397260;  % y1(t) = at^2
b = 2.739726027397260;  % y3(t) = -b(t-1)^2+1
c = -2.029426686960933; % y2(t) = ct^3 + dt^2 + gt + h
d = 3.044140030441400;
g = -0.015220700152207;
h = 0.000253678335870;


t1 = 0.05*tP;   
t2 = 0.95*tP;
tt1 = 0:dt:tP;
for i=1:length(tt1)
    t = tt1(i);
    if (t <t1) 							
        r1(i) = a*(t/tP)^2;                             %a*pow((t/tP),2);                
    elseif ( (t >=t1 ) && ( t<t2) )  
        r1(i) = c*(t/tP)^3 + d*(t/tP)^2 + g*(t/tP) + h; %c*pow((t/tP),3) + d*pow((t/tP),2) + g*(t/tP) + h;
    elseif (t>=t2) 
        r1(i) = -b*( (t/tP) - 1.0 )^2 + 1;                  %-b*pow( ((t/tP) - 1),2) + 1;
    end
end


tt2 = tP:dt:tP+2*offset;
tt4 = (2*tP+2*offset):dt:(3*tP+2*offset);
r2 = ones(1,length(tt2));
r4 = zeros(1,length(tt4));
%%RHS Expansion Phase 
tt3 = (tP+2*offset):dt:(2*tP+2*offset);
r3 = r1 - 0.5;
r3 = -r3 + 0.5;

%plot(tt1,r1,'r.','MarkerSize',6); hold on;
%plot(tt2,r2,'r.','MarkerSize',6); hold on;
%plot(tt3,r3,'r.','MarkerSize',6); hold on;
%plot(tt4,r4,'r.','MarkerSize',6); hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LHS
%
%%%%%%%%%%%%%%%%%%%%%%%%%

tprev = offset;
t1 = 0.05*tP+offset;   
t2 = 0.95*tP+offset;
frac= 0.9;
ttt1 = offset:dt:(tP+offset);
for i=1:length(ttt1)
    t = ttt1(i);
    if (t <t1) 							
        Lgeo(i) =  a*( (t-tprev) /tP)^2;                             %a*pow((t/tP),2);                
    elseif ( (t >=t1 ) && ( t<t2) )  
        Lgeo(i) = c*( (t-tprev)/tP)^3 + d*( (t-tprev)/tP)^2 + g*( (t-tprev)/tP) + h; %c*pow((t/tP),3) + d*pow((t/tP),2) + g*(t/tP) + h;
    elseif (t>=t2) 
        Lgeo(i) = -b*( ((t-tprev)/tP) - 1.0 )^2 + 1;                  %-b*pow( ((t/tP) - 1),2) + 1;
    end
end

tt0 = 0:dt:offset;
L0 = zeros(1,length(tt0));

%Contraction
L1 = frac*Lgeo;

%Phase Two -> Expansion
ttt2 = (tP+offset):dt:(2*tP+offset);
L2 = Lgeo - 0.5;
L2 = -1.2*L2 + 0.4 - (1-frac);

%Phase Three -> Contraction to equilibrium
ttt3 = (2*tP+offset):dt:(3*tP+offset);
L3 = Lgeo - 1.0;
L3 = (0.3)*L3;


YAXISy = -0.75:dt:1.6;
YAXISx = zeros(1,length(YAXISy));

XAXISx = -0.25:dt:1.5;
XAXISy = zeros(1,length(XAXISx));

plot(tt1,r1,'r.','MarkerSize',8); hold on;
plot(tt0,L0,'.','MarkerSize',8); hold on;

plot(XAXISx,XAXISy,'k-','LineWidth',5);
plot(YAXISx,YAXISy,'k-','LineWidth',5);

plot(tt1,r1,'r.','MarkerSize',8); hold on;
plot(tt0,L0,'.','MarkerSize',8); hold on;

plot(tt2,r2,'r.','MarkerSize',8); hold on;
plot(tt3,r3,'r.','MarkerSize',8); hold on;
plot(tt4,r4,'r.','MarkerSize',8); hold on;

plot(ttt1,L1,'.','MarkerSize',8); hold on;
plot(ttt2,L2,'.','MarkerSize',8); hold on;
plot(ttt3,L3,'.','MarkerSize',8); hold on;
axis([-offset 3*tP+3*offset -0.75 1.5]);

vL1y = -0.75:0.25:0;
vL1x = offset*ones(1,length(vL1y));

vL2y = -0.75:0.25:0.9;
vL2x = (tP+offset)*ones(1,length(vL2y));

vL3y = -0.75:0.25:-0.21;
vL3x = (2*tP+offset)*ones(1,length(vL3y));

vL4y = -0.75:0.25:0.0;
vL4x = (3*tP+offset)*ones(1,length(vL4y));

vL5y = -0.75:0.25:0.0;
vL5x = (3*tP+2*offset)*ones(1,length(vL4y));

vR1y = 1:0.2:1.6;
vR1x = tP*ones(1,length(vR1y)); 

vR2y = 1:0.2:1.6;
vR2x = (tP+2*offset)*ones(1,length(vR2y)); 

vR3y = 0:0.2:1.6;
vR3x = (2*tP+2*offset)*ones(1,length(vR3y));

plot(vR1x,vR1y,'k--','LineWidth',1); hold on;
plot(vR2x,vR2y,'k--','LineWidth',1); hold on;
plot(vR3x,vR3y,'k--','LineWidth',1); hold on;

plot(vL1x,vL1y,'k--','LineWidth',1); hold on;
plot(vL2x,vL2y,'k--','LineWidth',1); hold on;
plot(vL3x,vL3y,'k--','LineWidth',1); hold on;
plot(vL4x,vL4y,'k--','LineWidth',1); hold on;
plot(vL5x,vL5y,'k--','LineWidth',1); hold on;
axis([-2*offset 3*tP+4*offset -1.0 1.75]);


legend('RIGHT','LEFT');