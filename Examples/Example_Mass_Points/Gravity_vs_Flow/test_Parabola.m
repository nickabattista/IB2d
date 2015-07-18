function test_Parabola

Lx = 1;
w = 0.2;

y = 0.41:0.002:0.59;

yVals = ( (Lx/2+w/2) - y ).*( (Lx/2-w/2) - y );

plot(yVals,y,'*'); hold on;
axis square;