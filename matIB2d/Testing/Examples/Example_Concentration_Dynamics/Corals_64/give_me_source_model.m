function[fs]=give_me_source_model(CL,Nc)

% Model for photosynthesis
if Nc==1
    fs=0.0215*CL(:,2);
elseif Nc==2
    fs=-0.0215*CL(:,2);
else
    'OOPS'
    
end