function test_Script()

% GEOMETRIC PARAMETERS
pi = 4*atan(1);
L1 = 8;                              % length of computational domain (m)
N1 = 512;                            % number of Cartesian grid meshwidths at the finest level of the AMR grid
bell_length = 2;                     % bell length (m)
bell_circumference = pi;
npts_bell = ceil(2*(bell_length/L1)*N1);  % number of pos along the length of the bell
npts_circ = 1; %number of pos along the circumference (if in 3D)
npts = npts_bell*npts_circ;	    % total number pos
ds1 = bell_length/(npts_bell-1);   % mesh spacing(m) along length of bell

% Values from Alben, Peng, and Miller
betao = 0.5;
betam = 0.3;
to = 0.4;
Zs = L1/8;
xRef = L1/2;

t=0:0.025:5;

for i=1:length(t) 


    %These are used to keep track of cycle number and time into the cycle
    pulse_time = t(i)-floor(t(i)); % determine time since beginning of first pulse
		
    % GIVE BELL STATES
    if (pulse_time<to)  %contract bell
        beta = betao+(betam-betao)*(pulse_time/to);
    else                % expand bell
        beta = betam+(betao-betam)*((pulse_time-to)/(1-to));
    end
		
    %---------------------
    %Xb and Yb are calculated here and will be used to determine new curvatures
    s1=0;
    zl = Zs;
    ro = xRef;
    
    % Pre-allocate memory for speed
    Xb_lam = zeros(npts,1);
    Yb_lam = zeros(npts,1);
    
    %top of bell
    Xb_lam(1)=ro;
    Yb_lam(1)=zl;

    %right side of bell
    for s1 = 2:ceil(npts_bell/2)
        thetaj = -1.55*(1-exp(-(s1)*ds1/beta));
        zl = zl + ds1*sin(thetaj);
        ro = ro + ds1*cos(thetaj);
        Xb_lam(s1)=ro;
        Yb_lam(s1)=zl;
    end

    zl = Zs;
    ro = xRef;
    
    %left side of bell
    for s1 = (ceil(npts_bell/2))+1:npts_bell
        s2=s1-(ceil(npts_bell/2));
        thetaj = -1.55*(1-exp(-(s2)*ds1/beta));
        zl = zl + ds1*sin(thetaj);
        ro = ro - ds1*cos(thetaj);
        Xb_lam(s1)=ro;
        Yb_lam(s1)=zl;
    end
   
    plot(Xb_lam,Yb_lam,'.'); hold on;
    axis([0 8 0 8]);
    pause(0.001); 
    clf;
    
end