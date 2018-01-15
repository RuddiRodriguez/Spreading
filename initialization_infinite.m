function [np,MTarryocupationtemp ] = initialization_infinite(ocupationnumber,numberpb,MTarryocupationtempindex,MTarryocupation)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

     
    
    
    [~,colocunozero] = (find (ocupationnumber~=0));
    
    try 
     MTLindex = (ceil ((colocunozero(end))))-(0:(numberpb-0));   
%      MTLindexgpu=gpuArray  (MTLindex) ;
%      MTarryocupationtempindexgpu=gpuArray(MTarryocupationtempindex);
    C = sort(MTarryocupationtempindex(ismembc(MTarryocupationtempindex,sort(MTLindex(MTLindex>0)))));
    
    %C = tmp([~~diff(tmp),true]);
%     C = intersect(MTarryocupationtempindex,MTLindex(MTLindex>0));
%     CC = setdiff(MTarryocupationtempindex,MTLindex(MTLindex>0)); 
%     CCC= setdiff(1:max(C),[CC C]); 
    np = length (C); 
     
     
%     MTarryocupationtemp = MTarryocupation(MTarryocupation~=0);   
    catch err 
        ocupationnumber(1)=1; 
        [~,colocunozero] = (find (ocupationnumber~=0)); 
        MTLindex = (ceil ((colocunozero(end))))-(0:(numberpb-0)); 
        np = length (MTarryocupationtempindex(end-(numberpb-1):end)); 
    end 
    
    
    MTarryocupationtemp = nonzeros(MTarryocupation);              
    
    if length (MTarryocupationtemp)< 2
        
        np=np;
    end
   
end

