function [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  
    controldensityneww=(sum(MTarryocupation(MTarryocupation~=0))+sum(ocupationnumber(ocupationnumber~=0)))./new_pos ;
%controldensityneww=(bsxfun(@plus, MTocutemp, memocutemp))/new_pos ;

end