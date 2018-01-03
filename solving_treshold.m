clear
temp = [10:10:100 ];
for i=1:length(temp)
V =@(x,b,kon,cl,koff,scale,betat,np) b*((kon*cl)-koff*exp(((x*scale*betat)/np)));
t = 0:15;                       % s
koff = 1.5;                     % 1/s
kon = 0.16;                     % 1/um2s
cl = temp(i);                       % # molecules / um2
np=4;
b = 8./1000;                    % reaction zone in um 
scale= (1).*1e-9;               % Barrier
Temp = 1.381e-23.*307.15;       % Temperatue  kelvin
kappa = 5e-20;                  % Bending modulus
sigmai = 2e-7;                   % Membrane tension
betat = 1.0/Temp;
fun = @(x) V(x,b,kon,cl,koff,scale,betat,np);    % function of x alone
x (i) = fzero(fun,1e-19);
end
kappa = (x.^2)./(8.*pi.*sigmai)
figure ; plot (temp, x.*1e12)