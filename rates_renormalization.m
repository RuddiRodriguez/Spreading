function [koofre,rates] = rates_renormalization(sigma,kappa,Vm,r0_ini,R_ini,koof,ratesi,betat,np,sigmai)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%rate renomalization
F0=(2.*pi.*sqrt(2.*kappa.*sigma));%+2.*pi*1.63e9*16e-18*(Vm.*1e-6)*(log(R_ini/r0_ini ));                                      % F0 to pull a tube
scale= (2).*1e-9;                                                       % barrier height
%F0 = 2e-9; F0 test
koofre = koof.*exp(((F0.*scale.*betat)./np));                            % off rate renormalization
rates=ratesi;                                                           % New rates
if sigmai==0
    rates(1,1) =koof  ;
else
    rates(1,1) =koofre  ;
end
%     rates = rates;%./rates(1,2);                                            % Can be normalized rates to 1
end

