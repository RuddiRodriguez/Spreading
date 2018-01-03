function [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
  
    controldensityneww=(( nnz(MTarryocupation(2:end)))+sum(nonzeros(ocupationnumber(2:end))))./new_pos ;
%controldensityneww=(bsxfun(@plus, MTocutemp, memocutemp))/new_pos ;

end