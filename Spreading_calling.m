function [pos,times,arrayrates,MTarryocupation,ocupationnumber,vector,controldensity,vinterp,results, parameters]= Spreading_calling (ratesi,...
    kappa,sigmai,maxsimutime,npin,density,initubel,densityindex)
if nargin < 1 || isempty (kappa)
    
    ratesi = [1.5 25 10000 10000 1000 1000];
    %ratesi = [0.0050    0.0400   15.8167   15.8167    1.0000    1.0000];
    Pb = ratesi(1,2)./(ratesi(1,1)+ratesi(1,2));
    Pu = ratesi(1,1)./(ratesi(1,1)+ratesi(1,2));
%     DM = (Pb*1)+Pu*(0.8);
%     rateDM = DM/(0.008*0.008);
%     ratesi(1,5:6) =  rateDM;
    %ratesi = [1.3 295 20 20 ];%100 100 100];
    %     ratesi = [1.3 25 120 120]%100 100 100];
end

if nargin < 2 || isempty (kappa)
    
    kappa =[  5].*1e-20;
end

if nargin < 3 || isempty (sigmai)
    
    sigmai = 4e-6;
end

if nargin < 4
    
    maxsimutime =70;
end

if nargin < 5
    
    npin = 1;
end

if nargin < 6
    
    density = [ 4   ] ;                                                         %particles per um2
end
if nargin < 7
    
    initubel = 200 ;                                                         %particles per um2
end
if nargin < 8
    
    densityindex = [1] ;                                                         %particles per um2
end
if nargin < 9
    
    v = 0.05 ;                                                         %particles per um2
end
%% 


%% 



figure ;
numsi=10;
numsi=1;
results = cell (numsi,length(densityindex),length(kappa));
MTLent= cell (numsi,length(densityindex),length(kappa)); 
parameters = cell (numsi,length(densityindex),length(kappa)); 
vt= 0 ;
for j =1:length(kappa)
    j
    for k =1:length(densityindex)
       
        vinterp=0;
        for i=1:numsi
              % maxsimutime =80+(120-80)*rand(1,1);
            sigmai = 3e-7+(2e-6-3e-7)*rand(1,1);
            npin = randi([1 3],1,1);
            
            [pos,times,globalrate,arrayrates,MTarryocupation,ocupationnumber,vector,yy,controldensity,vinterp,R_ini,r0_ini] = membrane_position_MT_Infinit_family_reaction (ratesi,...
                kappa(j),sigmai,maxsimutime,npin,density(1),initubel,densityindex(k),v);
            
            i
           
            results {i,k,j} = vector;
            parameters {i,k,j} = [sigmai maxsimutime R_ini npin];
            parameters {i,k,j} = [sigmai maxsimutime R_ini npin r0_ini];
            assignin('base', 'results', results);
            assignin('base', 'MTLent',  MTLent);
            assignin('base', 'parameters', parameters)
        end
        %vt (vt<=0)=NaN;
%         vtmean(k,j)=nanmean(nanmean(vt));
%         stdmean (k,j)
    end
end
% [fitresult, gof] = createFit_power(vector(:,1),vector(:,2)-vector(1,2));
% figure;contourf(kappa,densityindex,vtmean);
% colorbar
%vt (vt<=0)=NaN;