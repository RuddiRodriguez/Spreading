function [vector,yy,vinterp] = sampling_data(times,pos,controldensity,inter,tmax)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
post = pos (1:inter:end);
T= times  (1:inter:end);
Den = controldensity(1:inter:end);
try
vector = [T' post Den];
catch err
    
    post = post (1:length(T));
    Den = Den (1:length(T));
    vector = [T' post Den];
    pos = pos(1:length(times));
end
timeinter = 0:0.5:tmax;
yy = interp1(times,pos(1:length(times)),timeinter);
vinterp = diff (yy)./diff(timeinter);
% yy=NaN;
% vinterp=NaN;

end

