function  figure_control(pos,count,times,controldensity,subplot1,subplot2,subplot3,subplot4)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

plot1=plot(times(1:count), pos(1:count,1), times(count),pos(count,1),'.r','Parent',subplot1);
hold on
plot2=plot(times(1:count), controldensity(1:count,1),...
    times(count),controldensity(count,1),'.r','Parent',subplot2);
%plot (times(count),v,'.g','Parent',subplot3);
% plot(times(1:count), MTL(1:count,1)-pos(1:count,1), times(count),MTL(count,1)-pos(count,1),'.b','Parent',subplot4);

drawnow;



end

