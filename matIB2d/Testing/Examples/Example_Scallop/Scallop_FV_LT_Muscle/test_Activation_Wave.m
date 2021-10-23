function test_Activation_Wave()

xLag = 1:0.025:3;


freq = 10;                      % frequency of traveling wave down tube
Ltube = xLag(end ) - xLag(1);   % length of heart tube
buff = 0.25;                    % percent-buff on each end of heart tube
L_AR = (1-2*buff)*Ltube;        % length of the actual activation region of tube
SqWidth = L_AR/10;              % width of the square traveling activation wave 
v = (L_AR-SqWidth) * freq;      % traveling wave velocity
k = 2*pi*freq/v;                % Wave-number for traveling wave


time = 0:0.001:1;

for j=1:length(time)
    
    t = time(j);
    t = rem(t,1/freq);                % gives remainder after "modular arithmetic" ("fmod" in C++)

    xL = xLag(1)+ (buff)*Ltube  + (v*t);
    xR = xLag(1)+ (buff)*Ltube  + SqWidth + (v*t);
    xM = (xL+xR)/2;
    c = SqWidth/1.5;
    

    for i=1:length(xLag)

        x = xLag(i);
        if ( ( x >= xL ) && ( x <= xR ) )
            af_Val(i) = 1;                % Traveling Square Wave
            %af_Val = exp( -(x-xM)^2 / (2*c^2) ); % Traveling Gaussian Wave
        else
            af_Val(i) = 0.0;
        end

    end
    
    plot(xLag,af_Val,'r-'); hold on;
    plot(xLag,af_Val,'*'); hold on;
    strTitle = strcat('Time(s) = ',num2str(t));
    title(strTitle);
    axis([0.8 3.2 -0.2 1.2]);
    pause(0.1);
    clf;

end