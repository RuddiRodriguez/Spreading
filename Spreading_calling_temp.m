function [pos,times,globalrate,arrayrates,MTarryocupation,ocupationnumber,vector,interpovar,controldensity,vinterp]= Spreading_calling (simun,...
    tensionin,kappaini,ratesini,maxsimutime,npin,density,densityindex,Initubeli,plotflag)

if nargin < 1
    
    simun = 1;
end

if nargin < 2
    tensionin = [ 1 ].*2e-7;
   
end
if nargin < 3
    
   kappaini = [5].*1e-20;
end

if nargin < 4
    
   ratesini = [1.4 12 100 100 100 100];
end


if nargin < 5
    
   maxsimutime = 32;
end
if nargin < 6
    
   npin = 4;
end
if nargin < 7
    
   density = 1;
end
if nargin < 8
    
   densityindex = 12;
end
if nargin < 9
    
   initubel = 200;
end
if nargin < 10
    
   plotflag = 1;
end
vt= 0 ;

for i=1:simun
[pos,times,globalrate,arrayrates,MTarryocupation,ocupationnumber,vector,interpovar,controldensity,vinterp] = membrane_position_MT_Infinit_family_reaction(ratesini,kappaini,...
                                                                                                                tensionin,maxsimutime,npin,density,densityindex,Initubeli,plotflag);
                                                                                                            
                                                                                                            
vt(1:length(vinterp),i) = vinterp;

end
vt (vt<=0)=NaN;
end