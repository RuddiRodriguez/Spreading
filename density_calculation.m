function [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    MTocutemp = MTarryocupation(2:end);
    memocutemp = ocupationnumber(2:end);
    controldensityneww=(sum(MTocutemp(MTocutemp~=0))+sum(memocutemp(memocutemp~=0)))/new_pos ;


