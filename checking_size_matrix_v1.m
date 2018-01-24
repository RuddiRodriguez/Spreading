function [ocupationnumbertemp,MTarryocupationtemp,arrayratestemp,iMTLsize] = checking_size_matrix_v1(ocupationnumber,MTarryocupation,arrayrates)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    
    ocupationnumbertemp = zeros (1,length(ocupationnumber));
    ocupationnumbertemp(1:5) = ocupationnumber(1,end-4:end);
    %ocupationnumbertemp(1,length(ocupationnumber)+1:end) = 0;
    
    MTarryocupationtemp = zeros (1,length( MTarryocupation));
     MTarryocupationtemp(1:5) =MTarryocupation(1,end-4:end);
    %ocupationnumbertemp(1,length(ocupationnumber)+1:end) = 0;
%     
    arrayratestemp = zeros (6,length(arrayrates));
    arrayratestemp(:,1:5) = arrayratestemp(:,end-4:end);
%     arrayratestemp(:,length(arrayrates)+1:end) = 0;
%     
%     ocupationnumbertemp=ocupationnumbertemp(1,1+lastnonzeromembranes(end)-50:end);
%     MTarryocupationtemp=MTarryocupationtemp(1,1+lastnonzeromembranes(end)-50:end);
%     arrayratestemp=arrayratestemp(:,1+lastnonzeromembranes(end)-50:end);
    
%     
    

%     ocupationnumbertemp=ocupationnumbertemp;
%     MTarryocupationtemp = MTarryocupationtemp;
%     arrayratestemp = arrayratestemp;


 MTarryocupationtempindex = find(MTarryocupation~=0);
    if isempty (MTarryocupationtempindex)
        pos (end) = 0 ;
        tess=1;
    end


iMTLsize=length(ocupationnumbertemp);
end

