function [globalrate] = globalrate_calculation_family_reaction(globalrate,arrayrates)
%UNTITLED9 Summary of this function goes here

%   Detailed explanation goes here
% arrayratesnozeros = arrayrates(:);
% if size (globalrate,2)==2
%     globalrate(1,1) = (sum(arrayrates(1,:)));                             %off rate
%     globalrate(1,2) = (sum(arrayrates(2,:)));                             %on-rate
%        
% end
% if size (globalrate,2)==3
%     globalrate(1,1) = (sum(arrayrates(1,:)));                             %off rate
%     globalrate(1,2) = (sum(arrayrates(2,:)));                             %on-rate
%     globalrate(1,3) = (sum(arrayrates(3,:)));                             %membrane_diffusion
%          
% end
% 
% if size (globalrate,2)==4
%     globalrate(1,1) = (sum(arrayrates(1,:)));                             %off rate
%     globalrate(1,2) = (sum(arrayrates(2,:)));                             %on-rate
%     globalrate(1,3) = (sum(arrayrates(3,:)));                             %membrane_diffusion
%     globalrate(1,4) = (sum(arrayrates(4,:)));                             %membrane_diffusion
%     
%          
% end
% if size (globalrate,2)==5
%     globalrate(1,1) = (sum(arrayrates(1,:)));                             %off rate
%     globalrate(1,2) = (sum(arrayrates(2,:)));                             %on-rate
%     globalrate(1,3) = (sum(arrayrates(3,:)));                             %membrane_diffusion
%     globalrate(1,4) = (sum(arrayrates(4,:)));                             %MTmembrane_diffusion
%     globalrate(1,5) = (sum(arrayrates(5,:)));
%     
%          
% end
% if size (globalrate,2)==6
%     globalrate(1,1) = (sum(arrayrates(1,:)));                             %off rate
%     globalrate(1,2) = (sum(arrayrates(2,:)));                             %on-rate
%     globalrate(1,3) = (sum(arrayrates(3,:)));                             %membrane_diffusion
%     globalrate(1,4) = (sum(arrayrates(4,:)));                             %membrane_diffusion
%     globalrate(1,5) = (sum(arrayrates(5,:)));                             %MTmembrane_diffusion
%     globalrate(1,6) = (sum(arrayrates(6,:)));                             %MTmembrane_diffusion
    globalrate = (sum(arrayrates,2))';  
    
         
% end
    
     
     
end

