%Initialization

t = 1;                       % s
koff = 1.5;                     % 1/s
kon = 0.16;                     % 1/um2s
cl = 100;                       % # molecules / um2
D = 2;                          % um2/s  
b = 8./1000;                    % reaction zone in um 
scale= (1).*1e-9;               % Barrier
Temp = 1.381e-23.*307.15;       % Temperatue  kelvin
kappa = 5e-20;                  % Bending modulus
sigmai = 2e-7;                   % Membrane tension
R_ini  = 10.*1e-6;% + (15e-6-5e-6).*rand(1,100);
r0_ini  = sqrt(kappa./(2.*sigmai));
Lc = (((Temp)./(4.*pi.*kappa)).*((R_ini.^2)./r0_ini)).*1e6;
np = 2;
pos = 0;
V= 0 ;
koofre=0;
Vt=0;


betat = 1.0/Temp;
beta_ini = ((4.*pi.*kappa).*betat).*(r0_ini./( R_ini.^2));          
sigma = sigmai.*exp (beta_ini.*(pos(end).*1e-6)); 
F0=(2.*pi.*sqrt(2.*kappa.*sigma));%+2.*pi*1.63e9*16e-18*(V(i).*1e-6)*(log(R_ini/r0_ini )); 
koofre = koff.*exp(((F0.*scale.*betat)./np(1))); 
%%
%Equation
V = b.*((kon.*cl(1))-koofre);
pos (i) = pos (i-1)+V(i).*t(i);
x=0:10

%Equation
p = (1/sqrt(4.*pi.*D.*t)).*exp(-(((x-(V.*t)).^2)/(4.*D.*t)))


