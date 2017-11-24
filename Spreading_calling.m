function [pos,times,globalrate,arrayrates,MTarryocupation,ocupationnumber,vector,interpovar,controldensity,vinterp,vt,vtmean,results,kappa,density,densityindex]= Spreading_calling (ratesi,...
    kappa,sigmai,maxsimutime,npin,density,initubel,densityindex)
if nargin < 1 || isempty (kappa)
    
    ratesi = [1.5 12 16200 16200 1000 1000];
    %ratesi = [1.3 295 20 20 ];%100 100 100];
    %     ratesi = [1.3 25 120 120]%100 100 100];
end

if nargin < 2 || isempty (kappa)
    
    kappa =[  1].*1e-20;
end

if nargin < 3 || isempty (sigmai)
    
    sigmai = 2e-7;
end

if nargin < 4
    
    maxsimutime =1;
end

if nargin < 5
    
    npin = 4;
end

if nargin < 6
    
    density = [  1   ] ;                                                         %particles per um2
end
if nargin < 7
    
    initubel = 200 ;                                                         %particles per um2
end
if nargin < 8
    
    densityindex = [10] ;                                                         %particles per um2
end

figure ;

results = cell (3,length(densityindex),length(kappa))
vt= 0 ;
for j =1:length(kappa)
    j
    for k =1:length(densityindex)
        k
        vinterp=0;
        for i=1:1
            [pos,times,globalrate,arrayrates,MTarryocupation,ocupationnumber,vector,interpovar,controldensity,vinterp] = membrane_position_MT_Infinit_family_reaction (ratesi,...
                kappa(j),sigmai,maxsimutime,npin,density(1),initubel,densityindex(k));
            
            i
            vt(1:length(vinterp),i) = vinterp;
            results {i,k,j} = vector;
            
        end
        vt (vt<=0)=NaN;
        vtmean(k,j)=nanmean(nanmean(vt));
    end
end
% figure;contourf(kappa,densityindex,vtmean);
% colorbar
%vt (vt<=0)=NaN;