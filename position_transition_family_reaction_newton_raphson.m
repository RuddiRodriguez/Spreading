%function [mnumber] = position_transition_family_reaction_newton_raphson (arrayrates,globalrate,r1)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
r1 = r1 (1)
globalrate = globalrate (1)
arrayratesnozeros = arrayrates(:);
 arrayratesnozeros = arrayratesnozeros(arrayratesnozeros~=0);

    
   % mnumber = 1;
  epsilon = 0.02;
  sumalpha=0;
low=1;
high =length(arrayratesnozeros);
 mnumber = fix(length(arrayratesnozeros)/2.0);
  
    
    
    %sumalpha = arrayratesnozeros(1)/globalrate;
    while abs((r1)-sumalpha) >= epsilon
       % mnumber = mnumber + 1;
      sumalpha=(sum(arrayratesnozeros(1:mnumber))/globalrate);
      sumalphatemp = (-r1+(2*sumalpha))/(-1);
        
%         
    end
%end

