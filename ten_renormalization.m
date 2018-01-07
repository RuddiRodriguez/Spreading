function [sigma] = ten_renormalization(sigmai,kappa,betat,r0_ini,R_ini,pos)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%% Tension renormalization
    if (sigmai==0) %|| (pos (end) < Lcritical)
        sigma=sigmai;
    else
        beta_ini = ((4.*pi.*kappa).*betat).*(r0_ini./( R_ini.^2));           % renormalization factor
        sigma = sigmai.*exp (beta_ini.*(pos(end).*1e-6));                   % tension renormalization
    end
    %       sigma=sigmai;
end

