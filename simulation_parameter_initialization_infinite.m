function [koof,R_ini,r0_ini,sigmao,koofreo,new_pos,times,t_final,count,npi,sigma,Lcritical,t,pos,betat] = simulation_parameter_initialization_infinite(ratesi,maxsimutime,numberpb,sigmai,kappa,Tubeli)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Temp = 1.381e-23.*307.15;                                                   % Simulation temperature
betat = 1.0/Temp;                                                            % Inverse temperature

koof = ratesi(1,1);                                                         % off-rate

R_ini  = 10e-6+((15e-6)-(3e-6))*rand(1,1);                                   % Vesicles radius

r0_ini  = sqrt(kappa./(2.*sigmai));                                         % Tubes radius

 pos = Tubeli/1000;                                                              % Initial position
sigmao = sigmai;
koofreo =koof ;

new_pos = pos(end);                                                         %Initialization membrane position

t = 0;                                                                      % time0

times = 0;

t_final=maxsimutime;                                                        %T final

count = 0 ;                                                                 %Parameter controlling the loop number

npi=numberpb;                                                               % number of petides at the leading edge without zippered

sigma = sigmai;                                                             % Tension initial

Lcritical = ((Temp./(4.*pi.*kappa))...                                      %Critical length to force changing with tube length
    .*((R_ini.^2)./r0_ini)).*1e6;

% rng(sum(1000000*clock),'v4');

%   subplot(1,2,2);

% plot(MTarryocupationplot,'s','MarkerSize',6);hold on
% %     xlabel('number of X molecules');
% %     ylabel('number of Y molecules');
% axis( [0 30 0 1.1]);
end

