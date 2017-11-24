function [MTarryocupation,ocupationnumber,arrayrates,TubeLi,iMTLsize,iarraysize] = simulation_initialization_matrix_infinite(reactionn,density,initubel,densityindex)
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



MTLi = 10000;                                                               %Microtbule length in nm

iMTLsize = MTLi/8;                                                          % Microtubule discretization

MTarryocupation = (zeros (1,iMTLsize));                                       % Ocupation number intitialization MT
% MTarryocupationplot = ones (1,iMTLsize);

                                          % Initial number of peptides attached to the MT

ocupationnumber = (zeros(1,iMTLsize));                                        % Initial number  of peptides in the reaction area on the tube

 indexocup = randi ([1 (iarraysize)],1,densityindex);
 ocupationnumber (indexocup) = randi ([0 density],1,length(indexocup));
 ocupationnumber (iarraysize+1) = randi ([1 density],1,1); 
 ocupationnumber (1) = 1;
 MTarryocupation (indexocup) = 1;
% MTarryocupation (1:iarraysize) = 1;
  MTarryocupation (1) = 1;
 
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

