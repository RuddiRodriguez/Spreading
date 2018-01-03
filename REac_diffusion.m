% %% 
% % Solution of the reaction diffusion equation
% % 
% % $$p\text{?}\left(x,t\right)=\frac{1}{\sqrt{4\pi \mathrm{Dt}}}e^{-\frac{x^2 
% % }{4\mathrm{Dt}}}$$
% 
% x = -10:0.1:10 ;
% t = 0:6;
% D = 0.4; 
% %%
% close all
% % Create figure
% figure1 = figure('Name','Figure');
% 
% % Create axes
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
% for i = 1:length(t)
% p = (sqrt(1/(4.*pi.*D.*t(i)))).*exp(-x.^2./(4.*D.*t(i)));
%  plot (x,p);hold on
% end
% %%
% 
% % Create ylabel
% ylabel('p(x,t)','FontName','Avenir','FontSize',18);
% 
% % Create xlabel
% xlabel('Position x','FontName','Avenir','FontSize',18);
% 
% box(axes1,'on');
% % Set the remaining axes properties
% set(axes1,'FontName','Avenir','FontSize',14);
% 
% Solution of the reaction diffusion equation with drift

% $$p\text{?}\left(x,t\right)=\frac{1}{\sqrt{4\pi \text{Dt}}}e^{-\frac{{\left(x-\mathrm{Vt}\right)}^2 
% }{4\text{Dt}}}$$
%
x = 0:0.01:5 ;
t = [1 5 10];
% t = 0:50;                       % s
koff = 1.5;                     % 1/s
kon = 0.16;                     % 1/um2s
cl = 120;                       % # molecules / um2
D = 0.03;                          % um2/s  
b = 8./1000;                    % reaction zone in um 
scale= (1).*1e-9;               % Barrier
Temp = 1.381e-23.*307.15;       % Temperatue  kelvin
kappa = 5e-20;                  % Bending modulus
sigmai = 2e-7;                   % Membrane tension
R_ini  = 10.*1e-6 + (15e-6-5e-6).*rand(1,100);
r0_ini  = sqrt(kappa./(2.*sigmai));
Lc = (((Temp)./(4.*pi.*kappa)).*((R_ini.^2)./r0_ini)).*1e6;
np = 4;
pos = 0;
V= 0 ;
koofre=0;
V =0;
%%
close all
% Create figure
figure2 = figure('Name','Figure');

% Create axes
axes2 = axes('Parent',figure2);
hold(axes2,'on');
for i = 1:length(t)
 betat = 1.0/Temp;
beta_ini = ((4.*pi.*kappa).*betat).*(r0_ini./( R_ini(1).^2));          
% sigma = sigmai.*exp (beta_ini.*(pos(end).*1e-6)); 
F0=(2.*pi.*sqrt(2.*kappa.*sigmai));%+2.*pi*1.63e9*16e-18*(V(i).*1e-6)*(log(R_ini/r0_ini )); 
koofre = koff.*exp(((F0.*scale.*betat)./np(1))); 
V(i) = b.*((kon.*cl)-koofre);
% pos (i) = pos (i-1)+V(i).*t(i);
pd = (sqrt((1/(4.*pi.*D.*t(i))))).*exp(-((x-(V(i).*t(i))).^2)./(4.*D.*t(i)));
 plot (x,pd);hold on
end
%%
% Create ylabel
ylabel('p(x,t)','FontName','Avenir','FontSize',18);

% Create xlabel
xlabel('Position x','FontName','Avenir','FontSize',18);

box(axes2,'on');
% Set the remaining axes properties
set(axes2,'FontName','Avenir','FontSize',14);