function [pos,times,globalrate,arrayrates,MTarryocupation,ocupationnumber,vector,interpovar,controldensity,vinterp] = membrane_position_MT_Infinit_family_reaction (ratesi,...
    kappa,sigmai,maxsimutime,npin,density,initubel,densityindex)

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
 
numberpb = npin;
Tubeli=initubel;
[koof,R_ini,r0_ini,~,~,new_pos,times,t_final,count,~,~,Lcritical,t,pos,betat] = simulation_parameter_initialization_infinite(ratesi,maxsimutime,numberpb,sigmai,kappa,Tubeli);
%post (1)= pos;
ouft=0;inft=0;ctr=0;transitionkindtest=0;positiontrantest=0;
positiontran=0;
tranflag=0;

%% Matrix Initialization
reactionn = 6;%length (ratesi(ratesi~=0));
globalrate = zeros (1,reactionn);
tau = zeros (1,reactionn);
sigmao=ratesi (1,1);
%post = zeros (1000,1);
% currentFolder = pwd;
% 
% 
% mkdir (currentFolder,'results');
% f = fullfile(currentFolder,'results');

m=0;
Vm = 0.008*ratesi(1,2);
%[r1] = random_generator;
[MTarryocupation,ocupationnumber,arrayrates,Tubeli,iMTLsize,iarraysize,densitylb,densitylu] = simulation_initialization_matrix_infinite(reactionn,density, initubel,densityindex,r0_ini,ratesi);
[matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1] = matrix_transition_reaction_family();

MTocutemp = MTarryocupation(2:end);
    memocutemp = ocupationnumber(2:end);
   controldensity=sum(MTocutemp(MTocutemp~=0))+sum(memocutemp(memocutemp~=0)) ;
   densityluu=controldensity;
AL =(Tubeli/1000);
controldensity=controldensity/AL;
ocupacontrol = ocupationnumber(2);
%% Loop
%  spmd
%for i=1:1000000
while t <= t_final
    %% Simulation initialization
    
    count = count+1;                                                        %loop number
    
    MTarryocupationtempindex = find(MTarryocupation~=0);
    if isempty (MTarryocupationtempindex)
        pos (end) = 0 ;
        break;
    end
    [np,MTarryocupationtemp ] = initialization_infinite(ocupationnumber,numberpb,MTarryocupationtempindex,MTarryocupation);
    
    %
    %
    %
    
    %% Tension renormalization
    
    
    
    if (sigmai==0) %|| (pos (end) < Lcritical)
        sigma=sigmai;
    else
        
        
        beta_ini = ((4.*pi.*kappa).*betat).*(r0_ini./( R_ini.^2));           % renormalization factor
        sigma = sigmai.*exp (beta_ini.*(pos(end).*1e-6));                   % tension renormalization
    end
    %       sigma=sigmai;
    %%
    
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
    
    
    %% Rates Matrix initialization arrayrates is the propensity of the reaction
    %     ocupationnumber = gather (ocupationnumber);
     if count ==1
     [arrayrates] = arrayrates_values_family_reaction(count,MTarryocupationtemp,numberpb,iMTLsize,MTarryocupation,ocupationnumber,rates,koof,arrayrates,positiontran,tranflag,densitylb,densitylu);
    
     end
    
    %
    %% Global rate calculation
    
    % arrayratesnozeros = arrayrates(:);
    [globalrate] = globalrate_calculation_family_reaction(globalrate,arrayrates);
    
    %     prob = arrayrates./globalrate;                                       % Total propensity' that anything happens
    
    r1= rand(1,reactionn+1);
    
    %     test=ceil(0 + ((globalrate)-0)*rand(1,1000));
    %     figure;histogram(test);
    %
    %% Set time for simulation and ending
    
    %     tau(count) = exprnd(1/globalrate);
    [tau] = tau_calculation_family_reaction(tau,globalrate,r1);
    [mintau,Imintau] = min (tau);
    
    
   
    
   
    
    t = t + mintau;                                                     % update time
    
    if t >= t_final
        
        t = t_final;
        
        pos = [pos; pos(end)];
       % post (count+1) = pos;
        
        times = [times, t];
        
        controldensitynew=sum(MTarryocupation(MTarryocupation~=0))+sum(ocupationnumber(ocupationnumber~=0)) ;
   AL = 2*pi*(r0_ini*1000000)*(pos(end));
    controldensitynew=controldensitynew/( AL);
   controldensity = [controldensity; controldensitynew];
        break;
    end
    
    
    
    %% Tube progression
    
    %initialization
    
    
    if Imintau==1
        
        [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:),globalrate(1),r1(reactionn+1));
    end
    
    if Imintau==2
        [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:),globalrate(2),r1(reactionn+1));
    end
    
    if Imintau==3
        [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:), globalrate(3),r1(reactionn+1));
    end
    if Imintau==4
        [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:),globalrate(4),r1(reactionn+1));
    end
    if Imintau==5
        [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:),globalrate(5),r1(reactionn+1));
    end
    if Imintau==6
        [mnumber] = position_transition_family_reaction (arrayrates(Imintau,:),globalrate(6),r1(reactionn+1));
    end
    
    %% taking position and kind of transition
    if mnumber==0
        positiontran = 0;                                                   % integer part
        
        transitionkind = 0;                                                 % remaining part
        
        [ ocupationnumber,MTarryocupation,new_pos,pos,iarraysize,status] = transitions_no_EB_infinite_family_reaction_backup_6reaction(transitionkind,positiontran,ocupationnumber,MTarryocupation,MTarryocupationtemp,new_pos,pos,iarraysize,matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1 );
        
    end
    
    if mnumber~=0 
        
        positiontran = mnumber;
        
        transitionkind = Imintau;
%         [ ocupationnumbertt,MTarryocupationtt,new_postt,postt,iarraysizett,statustt,controlt ] = transitions_no_EB_infinite_family_reaction(transitionkind,positiontran,ocupationnumber,MTarryocupation,MTarryocupationtemp,new_pos,pos,iarraysize,...
%                                                                             matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1);
% lastnonzeromembranest = find (ocupationnumbertt>0);
%  lastnonzeroMTt = find (MTarryocupationtt>0);
% if lastnonzeromembranest(end)<=lastnonzeroMTt(end)
%         testb =1;
%     end 
        [ ocupationnumber,MTarryocupation,new_pos,pos,iarraysize,status,tranflag,ouft,inft,transitionkindtest,positiontrantest,ctr ] = transitions_no_EB_infinite_family_reaction_backup_6reaction(transitionkind,positiontran,ocupationnumber,...
            MTarryocupation,MTarryocupationtemp,new_pos,pos,iarraysize,matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1,ouft,inft ,ctr,transitionkindtest,positiontrantest);
       if count ~=1
           MTocutemp = MTarryocupation(2:end);
    memocutemp = ocupationnumber(2:end);
   controldensitynew=sum(MTocutemp(MTocutemp~=0))+sum(memocutemp(memocutemp~=0)) ;
   AL = (pos(end));
   controldensitynew=controldensitynew/( AL);
   controldensity = [controldensity; controldensitynew];
  densitylu = sum(memocutemp)/(AL);
  densitylb = sum(MTocutemp)/(AL);
     [arrayrates] = arrayrates_values_family_reaction(count,MTarryocupationtemp,numberpb,iMTLsize,MTarryocupation,ocupationnumber,rates,koof,arrayrates,positiontran,tranflag,densitylb,densitylu);
    
    end 
        
    end
    
    
    %     figure ; plot (ocupationnumber);hold on ; plot (MTarryocupation);
    
    %% Checking to finish the loop
    if status == 0
              
        break;
    end
    
%     if isempty ( ocupationnumber)
%         break;
%     end
    
     pos = [pos;new_pos];                                                    % update vector with all positions
     
     
     
     sigmao = [  sigmao ; koofre];                                            % update vector with  ar to check
% post (count+1)= new_pos ;  
% if count~=1 && pos(end)>pos(end-1)
%     transitionkind
% end
    times = [times, t]; 
    % update vector with all times
    %     koofreo = [koofreo ; koofre];
     MTocutemp = MTarryocupation(2:end);
    memocutemp = ocupationnumber(2:end);
   controldensitynew=sum(MTocutemp(MTocutemp~=0))+sum(memocutemp(memocutemp~=0)) ;
   AL = (pos(end));
   controldensitynew=controldensitynew/( AL);
   controldensity = [controldensity; controldensitynew];
   ocupacontrol = [ocupacontrol;ocupationnumber(2)];
   if (controldensity(end))>40
       testr = 1;
   end
   if pos(end)<=0
         pos = [pos;new_pos];
          times = [times, t];
          controldensity = [controldensity; controldensitynew];
        break;
    end
   
    %% Plotting
    
     if (mod(m,100000) == 0)
           subplot(1,3,1);
              plot(times(1:count), pos(1:count,1), times(count),pos(count,1),'.r');
% % % % % %         
% % % % % % %        plot(MTarryocupation(MTarryocupation~=0),'LineWidth',12);
    subplot(1,3,2);
            plot(times(1:count), controldensity(1:count,1),...
            times(count),controldensity(count,1),'.r');
 subplot(1,3,3);
% % % % % % %   plot(times(1:count), sigmao(1:count,1), times(count),sigmao(count,1),'.r');
%  ocupp = ocupationnumber(ocupationnumber~=0);
%  ocupp (ocupp==0)= NaN;
%  plot (ocupp,'ro');
% plot (times(count),ocupacontrol,'.r');
%  plot(times(1:count), ocupacontrol(1:count,1), times(count),ocupacontrol(count,1),'.r');
         drawnow;
% % % % 
%  [vector,yy] = sampling_data(times,pos,200,1);
% % vector = [times',pos,controldensity];
   fprintf('%i\n',times(end)); 
%  save(sprintf('%s/retemp',f),'yy','-ascii');

    end
    
    m=m+1;
    
end
%  end
%% Checking times
%  [post] = getting_var_smpd(pos);
%  [taut] = getting_var_smpd(tau);


if status == 0
    td = 0;
    while pos(end) >0.001 
        td = td +1;
        new_pos = pos(end) * exp (-0.00005*td);
        pos = [pos;new_pos];    
        
    end
%  figure;plot(pos)   
    
end 


 [vector,interpovar,vinterp] = sampling_data(times,pos,controldensity,50000,t_final);
end