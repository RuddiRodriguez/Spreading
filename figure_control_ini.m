function  [subplot1,subplot2,subplot3,subplot4]=figure_control_ini(figure1)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
subplot1 = subplot(2,2,1,'Parent',figure1);
hold(subplot1,'on');
% Create ylabel
ylabel({'x ({\mu}m)'});

% Create xlabel
xlabel({'Time (s)'},'FontSize',11,'FontName','Arial');

% Create title
title({'Microtubule and membrane position'});

box(subplot1,'on');

subplot2= subplot(2,2,2,'Parent',figure1);
hold(subplot2,'on');
% Create ylabel
ylabel({'{\rho}(molecules/{\mu}m^2)'});

% Create xlabel
xlabel({'Time (s)'});

% Create title
title({'Density'});
box(subplot2,'on');

subplot3=subplot (2,2,3,'Parent',figure1);
hold(subplot3,'on');

% Create ylabel
ylabel({'V ({\mu}m/min)'});

% Create xlabel
xlabel({'Time (s)'});

% Create title
title({'Microtubule growth rate'});
box(subplot3,'on');

subplot4=subplot(2,2,4,'Parent',figure1);
hold(subplot4,'on');
% Create ylabel
ylabel({'{\xi}'});

% Create xlabel
xlabel({'Time  (s )'});

% Create title
title({'Position(Mt)-Position(membrane)'});

box(subplot4,'on');


end

