function [arrayindexnon] = finding_zeros(ocupationnumber)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
count=0;
arrayindexnon=0;
for i = 1:length(ocupationnumber)
    if ocupationnumber(i)~=0
        count=count+1;
        arrayindexnon(count)=i;
        
    end
end
end

