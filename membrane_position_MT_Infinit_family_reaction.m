   
function [pos,times,globalrate,arrayrates,MTarryocupationttotal,ocupationnumbertotal,vector,yy,controldensity,vinterp,R_ini,r0_ini] = membrane_position_MT_Infinit_family_reaction (ratesi,...
    kappa,sigmai,maxsimutime,npin,density,initubel,densityindex,vm)

%%
% this function try to simulate the membrane spreading using
% input
%ratesi- initial values of on-rate and off-rate for a single molecule in  1/s
%kappa-  bending modulus in J
%sigmai- intila value of mebrane tension in N/m .
%
%
%% Input checking
close all
if nargin < 1 || isempty (ratesi)
    
    ratesi = [1.2 12 5 5 0 0];
    %ratesi = [1.3 295 20 20 ];%100 100 100];
    %     ratesi = [1.3 25 120 120]%100 100 100];
end

if nargin < 2 || isempty (kappa)
    
    kappa =5e-20;
end

if nargin < 3 || isempty (sigmai)
    
    sigmai = 2e-7;
end

if nargin < 4
    
    maxsimutime =8;
end

if nargin < 5
    
    npin = 4;
end

if nargin < 6
    
    density = 1 ;                                                         %particles per um2
end
if nargin < 7
    
    initubel = 32 ;                                                         %particles per um2
end
if nargin < 8
    
    densityindex = 1 ;                                                         %particles per um2
end
%% Parameters Intialization
%Number of peptides a tthe tip
numberpb = npin;
%Initial length of the tube
Tubeli=initubel;
%Several parameters  Initialization
[koof,R_ini,r0_ini,~,~,new_pos,times,t_final,count,~,~,Lcritical,t,pos,betat] = simulation_parameter_initialization_infinite(ratesi,maxsimutime,numberpb,sigmai,kappa,Tubeli);

positiontran=0;
tranflag=0;
order_P = 0;
figure1 = figure;
[subplot1,subplot2,subplot3,subplot4]=figure_control_ini(figure1);
reactionn = 6;%length (ratesi(ratesi~=0));
globalrate = zeros (1,reactionn);
tau = zeros (1,reactionn);
m=0;
Vm = 0.008*ratesi(1,2);
tempp=0;
counttotal=0;
ocupationnumbertotal=0;
MTarryocupationttotal=0;
% vm=0.09
%sigmao=ratesi (1,1);
%post = zeros (1000,1);
%[r1] = random_generator;
%% Folder_To save  the results
% currentFolder = pwd;
%
%
% mkdir (currentFolder,'results');
% f = fullfile(currentFolder,'results');

%% Matrix Initialization

[MTarryocupation,ocupationnumber,arrayrates,Tubeli,iMTLsize,iarraysize,densitylb,densitylu] = simulation_initialization_matrix_infinite(reactionn,density, initubel,densityindex,r0_ini,ratesi);
[matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1] = matrix_transition_reaction_family();
%% Density Calculation

MTocutemp = MTarryocupation(2:end);
memocutemp = ocupationnumber(2:end);
controldensity=sum(MTocutemp(MTocutemp~=0))+sum(memocutemp(memocutemp~=0));
AL =(Tubeli/1000);
controldensity=controldensity/AL;
controldensitynew=0;

%% Loop
% spmd
while t <= t_final
    %% Simulation initialization
    
    count = count+1;                                                        %loop number
    
    MTarryocupationtempindex = find(MTarryocupation~=0);
    if isempty (MTarryocupationtempindex)
        pos (end) = 0 ;
        break;
    end
    [np,MTarryocupationtemp ] = initialization_infinite_new(ocupationnumber(:),numberpb,MTarryocupationtempindex(:),MTarryocupation(:));
    
    %% Tension renormalization
    [sigma] = ten_renormalization(sigmai,kappa,betat,r0_ini,R_ini,pos);
    %% Rate Renormalization
    [~,rates] = rates_renormalization(sigma,kappa,Vm,r0_ini,R_ini,koof,ratesi,betat,np,sigmai);
    
    %% Rates Matrix initialization arrayrates is the propensity of the reaction
    %     ocupationnumber = gather (ocupationnumber);
    if count ==1
        [arrayrates] = arrayrates_values_family_reaction(count,MTarryocupationtemp,numberpb,iMTLsize,MTarryocupation,ocupationnumber,rates,koof,arrayrates,positiontran,tranflag,densitylb,densitylu,vm);
        
    end
    
    %% Global rate calculation
    
    % arrayratesnozeros = arrayrates(:);
    [globalrate] = globalrate_calculation_family_reaction(globalrate,arrayrates);
    r1= rand(1,reactionn+1);
    
    %% Set time for simulation and ending_finding the Family to fire
    
    %     tau(count) = exprnd(1/globalrate);
    [tau] = tau_calculation_family_reaction(tau,globalrate,r1);
    [mintau,Imintau] = min (tau);
    t = t + mintau;                                                     % update time
    if t >= t_final
        
        t = t_final;
        %order_P = [order_P;order_P_new];
        pos = [pos; pos(end)];
        times = [times, t];
        
        controldensitynew=sum(MTarryocupation(MTarryocupation~=0))+sum(ocupationnumber(ocupationnumber~=0)) ;
        AL = 2*pi*(r0_ini*1000000)*(pos(end));
        controldensitynew=controldensitynew/( AL);
        controldensity = [controldensity; controldensitynew];
        break;
    end
    %% Reaction inside a Family fired
    
    [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:),globalrate(Imintau),r1(reactionn+1));
    
    %% Transition
    if mnumber~=0
        
        positiontran = mnumber;
        
        transitionkind = Imintau;
        %  [controldensitynet] =density_calculation(MTarryocupation(2:end),ocupationnumber(2:end),new_pos);
        [ ocupationnumber,MTarryocupation,new_pos,pos,iarraysize,status,tranflag,controldensitynew,lastnonzeromembranes ] = transitions(transitionkind,positiontran,ocupationnumber,...
            MTarryocupation,MTarryocupationtemp,new_pos,pos,iarraysize,matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1,controldensity(end) );
        [controldensitynew] =density_calculation(MTarryocupation(2:end),ocupationnumber(2:end),new_pos);
        %             if controldensitynew ~=(controldensitynet)
        %                 test=1;
        %             end
        MTarryocupationtempindex = find(MTarryocupation~=0);
        if isempty (MTarryocupationtempindex)
            pos (end) = 0 ;
            break;
        end
        if (lastnonzeromembranes(end)==length(ocupationnumber))% || lastnonzeromembranes(end)>50)
            %             [ocupationnumber,MTarryocupation,arrayrates,iMTLsize,tempp] = checking_size_matrix(ocupationnumber,MTarryocupation,arrayrates,lastnonzeromembranes,tempp);
            %[ocupationnumber,MTarryocupation,arrayrates,iMTLsize] = checking_size_matrix_v1(ocupationnumber,MTarryocupation,arrayrates);
            %[ocupationnumber,MTarryocupation,arrayrates,iMTLsize,ocupationnumbertotal,MTarryocupationttotal] = checking_size_matrix_v1(ocupationnumber,MTarryocupation,arrayrates,ocupationnumbertotal,MTarryocupationttotal);
            %[ocupationnumbertemp,MTarryocupationtemp,arrayratestemp,iMTLsize,tempp] = checking_size_matrix(ocupationnumber,MTarryocupation,arrayrates,lastnonzeromembranes,tempp);
            counttotal=counttotal+1;
            [ocupationnumber,MTarryocupation,arrayrates,iMTLsize,ocupationnumbertotal,MTarryocupationttotal] = checking_size_matrix_v1(ocupationnumber,MTarryocupation,arrayrates,counttotal,ocupationnumbertotal,MTarryocupationttotal);
        end
        if count ~=1
            
            [arrayrates] = arrayrates_values_family_reaction(count,MTarryocupationtemp,numberpb,iMTLsize,MTarryocupation,ocupationnumber,rates,koof,arrayrates,positiontran,tranflag,densitylb,densitylu,vm);
            
        end
        
    end
    
    %% Checking to finish the loop
    if status == 5
        
        pos = [pos;new_pos];
        times = [times, t];
        controldensity = [controldensity; controldensitynew];
        
        break;
    end
    
    if isempty ( ocupationnumber)
        break;
    end
    
    if new_pos<=0
        pos = [pos;new_pos];                                                    % update vector with all positions
        times = [times, t];
        controldensity = [controldensity; controldensitynew];
        break;
    end
    
    
    %% Updating Outputs
    
    
    pos = [pos;new_pos];                                                    % update vector with all positions
    times = [times, t];
    controldensity = [controldensity; controldensitynew];
    %     if pos(end)>0.1
    %         ss=1;
    %     end
    
    
    
    %     if pos(end)<0
    %         new_pos = 0;
    %         pos = [pos;new_pos];
    %         order_P = [order_P;order_P_new];
    %         MTL = [MTL;new_posMTL,];
    %         times = [times, t];
    %         controldensity = [controldensity; controldensitynew];
    %
    % %         break;
    %
    %     end
    %
    %% Plotting
    
    if (mod(m,10000) == 0)
           figure_control(pos,count,times,controldensity,subplot1,subplot2,subplot3,subplot4);
        % figure (1);plot(pos(count),rates(1,1),'o');hold off
        %           drawnow;
        
        fprintf('%4.2f\n',times(end));
        fprintf('%4.2f\n',t_final-times(end));
        assignin('base', 'times', times);
        assignin('base', 'pos', pos);
        assignin('base', 'controldensity', controldensity);
        
    end
    
    m=m+1;
    
end
%   end
%% Checking times
%  [post] = getting_var_smpd(pos);
%  [taut] = getting_var_smpd(tau);


% if status == 5
%     td = 0;
%     while pos(end) >0.001
%         td = td +1;
%         new_pos = pos(end) * exp (-0.00005*td);
%         pos = [pos;new_pos];
%
%     end
%     %  figure;plot(pos)
%
% end


[vector,yy,vinterp,~] = sampling_data(times,pos,2,t_final);
end