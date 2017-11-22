function [tau] = tau_calculation_family_reaction(tau,globalrate,r1)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
% if  size (tau,2) ==2
%     
%     tau(1,1)=(1/globalrate(1,1)).*log(1/r1(1));
%     tau(1,2)=(1/globalrate(1,2)).*log(1/r1(2));
% end
% if  size (tau,2) ==3
%     tau(1,1)=(1/globalrate(1,1)).*log(1/r1(1));
%     tau(1,2)=(1/globalrate(1,2)).*log(1/r1(2));
%     tau(1,3)=(1/globalrate(1,3)).*log(1/r1(3));
% end
% if size (tau,2) ==4
%     tau(1,1)=(1/globalrate(1,1)).*log(1/r1(1));
%     tau(1,2)=(1/globalrate(1,2)).*log(1/r1(2));
%     tau(1,3)=(1/globalrate(1,3)).*log(1/r1(3));
%     tau(1,4)=(1/globalrate(1,4)).*log(1/r1(4));
% end
% if size (tau,2) ==5
%     tau(1,1)=(1/globalrate(1,1)).*log(1/r1(1));
%     tau(1,2)=(1/globalrate(1,2)).*log(1/r1(2));
%     tau(1,3)=(1/globalrate(1,3)).*log(1/r1(3));
%     tau(1,4)=(1/globalrate(1,4)).*log(1/r1(4));
%     tau(1,5)=(1/globalrate(1,5)).*log(1/r1(5));
% end
% 
% if size (tau,2) ==6
% %     tau(1,1)=(1/globalrate(1,1)).*log(1/r1(1));
% %     tau(1,2)=(1/globalrate(1,2)).*log(1/r1(2));
% %     tau(1,3)=(1/globalrate(1,3)).*log(1/r1(3));
% %     tau(1,4)=(1/globalrate(1,4)).*log(1/r1(4));
% %     tau(1,5)=(1/globalrate(1,5)).*log(1/r1(5));
% %     tau(1,6)=(1/globalrate(1,6)).*log(1/r1(6));
    tau=(1./globalrate).*log(1./r1(1:end-1));
% end
end

