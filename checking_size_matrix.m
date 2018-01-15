function [ocupationnumbertemp,MTarryocupationtemp,arrayratestemp,iMTLsize,tempp] = checking_size_matrix(ocupationnumber,MTarryocupation,arrayrates,lastnonzeromembranes,tempp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if lastnonzeromembranes(end)==length(ocupationnumber)-5
    ocupationnumbertemp = zeros (1,length(ocupationnumber)+20);
    ocupationnumbertemp(1,1:length(ocupationnumber)) = ocupationnumber;
    ocupationnumbertemp(1,length(ocupationnumber)+1:end) = 0;
    
    MTarryocupationtemp = zeros (1,length(MTarryocupation)+20);
    MTarryocupationtemp(1,1:length(MTarryocupation)) = MTarryocupation;
    MTarryocupationtemp(1,length(MTarryocupation)+1:end) = 0;
%     
    arrayratestemp = zeros (6,length(arrayrates)+20);
    arrayratestemp(:,1:length(arrayrates)) = arrayrates;
    arrayratestemp(:,length(arrayrates)+1:end) = 0;
else
    ocupationnumbertemp=ocupationnumber;
    MTarryocupationtemp = MTarryocupation;
    arrayratestemp = arrayrates;
end
if lastnonzeromembranes(end)>50 && lastnonzeromembranes(end)>tempp
    
    ocupationnumbertemp = zeros (1,length(ocupationnumber));
    ocupationnumbertemp(1,1:end-1) = ocupationnumber(1,2:end);
    %ocupationnumbertemp(1,length(ocupationnumber)+1:end) = 0;
    
    MTarryocupationtemp = zeros (1,length( MTarryocupation));
     MTarryocupationtemp(1,1:end-1) =MTarryocupation(1,2:end);
    %ocupationnumbertemp(1,length(ocupationnumber)+1:end) = 0;
%     
    arrayratestemp = zeros (6,length(arrayrates));
    arrayratestemp(:,1:end-1) = arrayratestemp(:,2:end);
%     arrayratestemp(:,length(arrayrates)+1:end) = 0;
%     
%     ocupationnumbertemp=ocupationnumbertemp(1,1+lastnonzeromembranes(end)-50:end);
%     MTarryocupationtemp=MTarryocupationtemp(1,1+lastnonzeromembranes(end)-50:end);
%     arrayratestemp=arrayratestemp(:,1+lastnonzeromembranes(end)-50:end);
    
%     
    
else
    ocupationnumbertemp=ocupationnumbertemp;
    MTarryocupationtemp = MTarryocupationtemp;
    arrayratestemp = arrayratestemp;
end



tempp=lastnonzeromembranes(end);
iMTLsize=length(ocupationnumbertemp);
end

