function test_Concentration()

X_1 = 0:0.025:1;
X_0 = 0:0.025:1;
for i=1:length(X_0)
    for j=1:length(X_1)
        Z(i,j) = 0.5*(tanh(80.0*(X_1(j)-0.5*X_1(end)-0.001*cos(2*pi*X_0(i))))+1);
    end
end


surf(X_0,X_1,Z)