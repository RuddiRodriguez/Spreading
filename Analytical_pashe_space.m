% Parameter initialization 
clear
t = 0:0.5:50;                       % s
koff = 1.5;                     % 1/s
kon = 0.16;                     % 1/um2s
cl = [90];                       % # molecules / um2
D = 2;                          % um2/s  
b = 8./1000;                    % reaction zone in um 
scale= (1).*1e-9;               % Barrier
Temp = 1.381e-23.*307.15;       % Temperatue  kelvin
kappa = [5].*1e-20;                  % Bending modulus
sigmai = 2e-7;                   % Membrane tension
R_ini  = 10.*1e-6 + (15e-6-5e-6).*rand(1,100);
r0_ini  = sqrt(kappa./(2.*sigmai));
Lc = (((Temp)./(4.*pi.*kappa)).*((R_ini.^2)./r0_ini)).*1e6;
np = 1;
pos = 0;
V= 0 ;
koofre=0;
Vt=0;
post=0;
%%
%Calculations
% figure(1) ;
count=0;
for j =1:length(R_ini)
    count = count+1;
    V=0;
    pos=0;
    counti=0;
for i=2:length (t)
    counti=counti+1;
betat = 1.0/Temp;
beta_ini = ((4.*pi.*kappa(1)).*betat).*(r0_ini(1)./( R_ini(j).^2));          
sigma = sigmai.*exp (beta_ini.*(pos(end).*1e-6)); 
F0=(2.*pi.*sqrt(2.*kappa(1).*sigma));%+2.*pi*1.63e9*16e-18*(V(i).*1e-6)*(log(R_ini/r0_ini )); 
koofre = koff.*exp(((F0.*scale.*betat)./np(1))); 
%%
%Equation
V(i) = b.*((kon.*cl(1))-koofre);
pos (i) = pos (i-1)+V(i).*t(i);
% if pos (i)<=0
%     pos (i)=0;
% end

end
% subplot (1,2,1)
%   plot (t(2:end),V(2:end));hold on;
%   subplot(1,2,2)
%   plot (t(2:end),pos(2:end));hold on;
 Vt (1:length(V),count)=V;
 post (1:length(pos),count)=pos;
 drawnow;
 
end
%% 
% for k = 1:size (Vt,2)
%    if Vt(2,k)>0
%        Is = Vt (2:end,k);
%      stgsv = [max(Is), 20, 3, min(Is)];
%     [estimates,eval] = fit_survival(t(2:end)', Is, stgsv);% MT reference
%     figure(1);plot(t(2:end),Is,'b-o',t(2:end),eval,'r-');hold on
% xlim([0 100])
% ylim([0 0.2])
%     
%     reftipr= estimates(2);
%     factor=abs(t-reftipr);
%     [col,ix(k)]=min (factor);
%    else
%        ix(k)=NaN;
%    end
%     
% end


%%
%plotting 
meanv=mean (Vt,2);
stdv = std (Vt,[],2)
 %errorbar(Exp1(:,1),Exp1(:,2),Exp1(:,3),'o');hold on
  errorbar(t(2:end),meanv(2:end),stdv(2:end));
 % plot(t(2:end),meanv(2:end));
%errorbar(Exp200(:,1),Exp200(:,2),Exp200(:,3),'o');
% xlim([0 60])
% ylim([0 0.2])

% A = b*kon;
% B = b*koff;
%  f = @ (cl,y) b*((kon*cl)-(koff.*exp(y)));
% % x0 = [500,1]; 
% % [xmin, fval] = fminsearch(@(x)f(x(1),x(2)), x0)
% %x = linspace(1,10,20);
% y = linspace(0,4,length(cl));
%  figure (2)
%  [CL,Y] = meshgrid(cl,y);
%  surf(CL,Y,f(CL,Y), 'EdgeColor','none', 'FaceAlpha',0.5)
% hold on
% plot3(xmin(1),xmin(2), f(xmin(1),xmin(2)), 'gp', 'MarkerSize', 10, 'MarkerFaceColor','g')
% hold off
% xlabel('X')
% ylabel('Y')

