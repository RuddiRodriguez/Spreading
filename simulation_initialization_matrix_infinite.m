function [MTarryocupation,ocupationnumber,arrayrates,TubeLi,iMTLsize,iarraysize,densitylb,densitylu] = simulation_initialization_matrix_infinite(reactionn,density,initubel,densityindex,r0_ini,ratesi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
                                                          % initial number of bounf peptides

TubeLi = initubel;                                                               % Intial tube length in nm by default 32 in nm

% pwidth = (25./13)./1000                                                     % width of on dimer um
% 
% Areatoi = pwidth.*(TubeLi./1000);                                           % area total um2
% 
% maxnumberpep = ceil(Areatoi.*densityi);                                      % max number of petides allowed 

iarraysize = ceil(TubeLi/8);                                                      % Initial discretization rates arrate



MTLi =400;                                                               %Microtbule length in nm

iMTLsize = (MTLi/8);                                                          % Microtubule discretization

MTarryocupation = zeros (1,iMTLsize);                                       % Ocupation number intitialization MT
% MTarryocupationplot = ones (1,iMTLsize);

                                          % Initial number of peptides attached to the MT
  densitys = 10000;


densityl = 2*pi*(r0_ini*1000000)*densitys;

Pb = ratesi(1,2)./(ratesi(1,1)+ratesi(1,2));
Pu = ratesi(1,1)./(ratesi(1,1)+ratesi(1,2));
densitylb = densityl*(Pb);
densitylu =densityl*(Pu);

numberb = (ceil(densitylb*(initubel/1000)))-0;
numberu = (ceil(densitylu*(initubel/1000)))-0;

 [numberbtemp] = random_fix_sum(numberb,iarraysize,1 );
 [numberutemp] = random_fix_sum(numberu,iarraysize,3 );
 


                                          
                                          
                                          
                                          

ocupationnumber = zeros(1,iMTLsize);                                        % Initial number  of peptides in the reaction area on the tube

 indexocup = randi ([2 (iarraysize)],1,length(numberutemp));
 if numberu>(iarraysize) %|| isnan(numberu)||numberu==0
 ocupationnumber (1:iarraysize) = 1;
 else
    ocupationnumber (indexocup) = numberutemp;%randi ([0 density],1,length(indexocup)); 
 end
 
 ocupationnumber (iarraysize+1) = randi ([3 3],1,1); 
 ocupationnumber (1) = 1;
 
 indexocupMT = randi ([1 (iarraysize)],1,length(numberbtemp));
 
 if numberb>(iarraysize)
      MTarryocupation (1:iarraysize) = 1;
 else
 
 MTarryocupation (indexocupMT) = 1;
 MTarryocupation (iarraysize) = 1;
  MTarryocupation (1) = 1;
 end

densitylb = sum( MTarryocupation(2:end))/(initubel/1000);
densitylu =sum(ocupationnumber(2:end) )/(initubel/1000);
%  
 
% ocupationnumber (1:iarraysize+1) = randi ([0 1],1,iarraysize+1);  
%%ocupationnumber (iarraysize) = randi ([1 3],1,iarraysize);                  % Initial number  of peptides in the reaction area on the tube
% maximuforindexarray = maxnumberpep-sum(ocupationnumber);
% if maximuforindexarray>iarraysize
%    maximuforindexarray=iarraysize; 
% end
% indextemp=randi ([0 1],1,maximuforindexarray);
% [colindex]=randperm(maximuforindexarray);
% MTarryocupation (colindex) = 1;     
% % sumocu=0;
% % count = 0 ;
% % while sumocu<maxnumberpep 
% %     count = count+1;
% %     MTarryocupation (count) = randi ([0 1],1,1); 
% %     sumocutemp = sum (MTarryocupation);
% %     sumocu = sumocutemp+sum(ocupationnumber);
% % end

% arrayrates = zeros (2,iMTLsize);                                            % Array of rates
arrayrates =(zeros (reactionn,iMTLsize));       
end

