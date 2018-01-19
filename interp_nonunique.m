function [your_vals] = interp_nonunique(x1,y1,x2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%[C,ia,idx] = unique(x1,'stable');
[C,ia,idx] = UniqueStable(x1);
val = accumarray(idx,y1,[],@mean); %You could take something other than the mean.
your_vals = interp1(C,val,x2,'linear','extrap'); %see interp1() for other interpolation methods. Extrapolation is dangerous.

end

