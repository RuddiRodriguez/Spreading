function [ocupationnumbertemp,MTarryocupationtemp,arrayratestemp,iMTLsize] = checking_size_matrix(ocupationnumber,MTarryocupation,arrayrates,lastnonzeromembranes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if lastnonzeromembranes(end)==length(ocupationnumber)-10
    ocupationnumbertemp = zeros (1,length(ocupationnumber)+50);
    ocupationnumbertemp(1,1:length(ocupationnumber)) = ocupationnumber;
    ocupationnumbertemp(1,length(ocupationnumber)+1:end) = 0;
    
    MTarryocupationtemp = zeros (1,length(MTarryocupation)+50);
    MTarryocupationtemp(1,1:length(MTarryocupation)) = MTarryocupation;
    MTarryocupationtemp(1,length(MTarryocupation)+1:end) = 0;
%     
    arrayratestemp = zeros (6,length(arrayrates)+50);
    arrayratestemp(:,1:length(arrayrates)) = arrayrates;
    arrayratestemp(:,length(arrayrates)+1:end) = 0;
else
    ocupationnumbertemp=ocupationnumber;
    MTarryocupationtemp = MTarryocupation;
    arrayratestemp = arrayrates;
end
iMTLsize=length(ocupationnumbertemp);
end

