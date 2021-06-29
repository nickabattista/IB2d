function test_HyperTangents()


xVec = 0:0.001:0.5;
sF=1;
coeff=40;

for i=1:length(xVec)
    
    x = xVec(i);
    
    if x < 0.1
        
        uX_Tar(i) = sF*( tanh(coeff*(0.05-x) ) ); 
    
    elseif x <= 0.2

        uX_Tar(i) = sF*( tanh(coeff*(x-0.15) ) ); 
        
    elseif x < 0.3
        

        uX_Tar(i) = sF*( tanh(coeff*(0.25-x) ) ); 

        
    elseif x<0.40
        
        uX_Tar(i) = sF*( tanh(coeff*(x-0.35) ) ); 

    else
        uX_Tar(i) = sF*( tanh(coeff*(0.45-x) ) ); 
    end

end

plot(xVec,uX_Tar,'.'); hold on;
