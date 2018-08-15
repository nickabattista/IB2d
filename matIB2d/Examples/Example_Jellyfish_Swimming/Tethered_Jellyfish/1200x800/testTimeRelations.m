function testTimeRelations()

current_time = 0:0.005:1.0;

for i=1:length(current_time)

period = 0.2;                                      %Period for 1 Cycle
t(i) = mod( current_time(i) , period);                     %Gives remainder time from full cycle to use for phases
percentOFFSET = 0.25;                              %Desired percentOFFSET for LEFT ARM
offset = percentOFFSET*period;                     %Gives offset for LEFT ARM
%tt(i) = ((period-offset)/offset)*t(i) - (period-offset); %TIME FOR LEFT ARM

%tt(i) = t(i) + (period-offset);

if t(i) < offset
    tt(i) = t(i) + (period-offset);
elseif ( ( t(i) >= offset ) && ( t(i) <= period ) )
    tt(i) = t(i) - offset;
end

end

plot(current_time,t,'*'); hold on;
plot(current_time,tt,'r*');
axis([-0.1 1.1 -0.1 0.3]);