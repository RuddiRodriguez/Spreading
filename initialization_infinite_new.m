function [np,MTarryocupationtemp ] = initialization_infinite_new(ocupationnumber,numberpb,MTarryocupationtempindex,MTarryocupation)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
MTarryocupationtemp = nonzeros(MTarryocupation);
if length (MTarryocupationtemp)< numberpb
   np = length (MTarryocupationtemp);
else
allIdx = 1:numel(MTarryocupation) ;
colocunozero = allIdx(MTarryocupation ~=0 ) ;
try
np=length(nonzeros(MTarryocupation(colocunozero(end)-(numberpb-1):colocunozero(end))));
catch err
    return;
end
end
% colocunozero = myFind_mex(ocupationnumber);
% MTLindex = (ceil ((colocunozero(end-1))))-(0:(numberpb-1));
% 
% C = (MTarryocupationtempindex(ismembc(MTarryocupationtempindex,sort(MTLindex(MTLindex>0)))));
% 
% np = length (C);
% if np<2
%    test=1
% end
% ocupationnumber(1)=1;
% [~,colocunozero] = (find (ocupationnumber~=0));
% MTLindex = (ceil ((colocunozero(end))))-(0:(numberpb-0));
% np = length (MTarryocupationtempindex(end-(numberpb-1):end));





end

