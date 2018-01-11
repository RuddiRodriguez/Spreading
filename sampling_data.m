function [vector,yy,vinterp] = sampling_data(times,pos,controldensity,inter,tmax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
post = pos (1:inter:end);
T= times  (1:inter:end);
Den = controldensity(1:inter:end);
try
vector = [T post Den];
catch err
    
    post = post (1:length(T));
    Den = Den (1:length(T));
    vector = [T' post Den];
    pos = pos(1:length(times));
end
timeinter = 0:3:tmax;
try
% yy = interp1(times,pos(1:length(times)),timeinter);
[yy,~,~,~]=bin_data_myy(times,pos(1:length(times)),timeinter);
catch err
% ve = cumsum(ones(size(pos))).*pos*eps                       % Scaled Offset For Non-Zero Elements
% ve = ve + cumsum(ones(size(pos))).*(pos==0)*eps             % Add Scaled Offset For Zero Elements
% vi = pos + ve  
%     pos = consolidator(pos);
[yy,~,~,~]=bin_data_myy(times,pos(1:length(times)),timeinter)
%      [~,indexuni]=unique(pos) ;
%     yy = interp1(times,pos(1:length(times)),timeinter);
% yy=NaN;
% return;
end
figure;plot (times,pos);hold on;plot(timeinter,yy,'ro');hold on
vinterp = diff (yy)./diff(timeinter)';
% yy=NaN;
% vinterp=NaN;

end

