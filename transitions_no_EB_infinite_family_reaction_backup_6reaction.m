function [ ocupationnumber,MTarryocupation,new_pos,pos,iarraysize,status,tranflag] = transitions_no_EB_infinite_family_reaction_backup_6reaction(transitionkind,positiontran,ocupationnumber,MTarryocupation,MTarryocupationtemp,new_pos,pos,iarraysize,...
    matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1)
%% Conditions for tubes growth , bound and unbound of peptides
%1- the position of the transition is given by the integer part  of the
% mnumber/number of transitions ---(2), the  transition is given by the remainig of division (remaiinig 0 bound, remaining 1 , unbound)
% in each position the first rate is bound the second rate is unbound
%(remaiinig 0 bound, remaining 1 , unbound)

%% Lattice



tranflag=0;
status=1;
lastnonzeromembranes = find (ocupationnumber>0);
lastnonzeroMT = find (MTarryocupation>0);
% if transitionkind == 0 &&  positiontran ==0
%
%     new_pos = pos(end) + 0;
%     iarraysize = iarraysize+0 ;
%     ocupationnumber  = ocupationnumber ;
%     MTarryocupation  = MTarryocupation ;
%
% end
% if transitionkind == 3
%
%     testt = 1;
% end


%bound at the lattice
if transitionkind == 2 &&     positiontran~=0 && positiontran < lastnonzeromembranes(end)-1 && positiontran ~=0 && MTarryocupation (positiontran)~=1
    new_pos = pos(end) + 0;
    iarraysize = iarraysize+0 ;
    %     ocupationnumber (positiontran) = 1;
    %     ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
    if positiontran~=1
        ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (2,:);
        MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (2,:);
    else
        ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (2,:);
        MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (2,:);
    end
    %     if ocupationnumber (positiontran)==0
    %         ocupationnumber (positiontran)=randi ([1 3],1,1);
    %     end
    % MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (2,:);
    % MTarryocupation (positiontran) = MTarryocupation (positiontran)+1;
    
    %  if lastnonzeromembranes(end)<lastnonzeroMT(end)
    %         testb =1
    %     end
    
end



%unbound at the lattice
try
    if  transitionkind == 1 &&    positiontran~=0 && positiontran < lastnonzeromembranes(end)-1  && MTarryocupation (positiontran)==1
        
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        %     ocupationnumber (positiontran) = ocupationnumber (positiontran)+1;
        %     MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
        if positiontran~=1
            ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (1,:);
            MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (1,:);
        else
            ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (1,:);
            MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (1,:);
        end
        
        %    if lastnonzeromembranes(end)<lastnonzeroMT(end)
        %         testb =1
        %     end
    end
catch err
    if isempty(lastnonzeromembranes)==1
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber  = ocupationnumber ;
        MTarryocupation  = MTarryocupation ;
    end
    
end


%% Edge


%bound at the tip edge
if transitionkind == 2 &&   positiontran~=0 && positiontran == lastnonzeromembranes(end) && MTarryocupation (positiontran)==0
    
    new_pos = pos(end) + 0.008;
    iarraysize = iarraysize+1 ;
    ocupationnumber(2:positiontran+1)=ocupationnumber(1:positiontran);
    %     ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
    %     ocupationnumber (positiontran+1) = randi ([1 1],1,1);
    %     if ocupationnumber (positiontran)==0
    %         ocupationnumber (positiontran)=randi ([1 3],1,1);
    %     end
    
    MTarryocupation (positiontran) = MTarryocupation (positiontran)+1;
    
    
    %   if lastnonzeromembranes(end)<lastnonzeroMT(end)
    %         testb =1
    %     end
    
end

% if positiontran == length ((ocupationnumber(ocupationnumber~=0)))-1 &&transitionkind == 1 && MTarryocupation (positiontran-1)==1
%     return ;
% end
try
    %unbound at the tip edge
    if transitionkind == 1 &&  positiontran == lastnonzeromembranes(end)-1 && MTarryocupation (positiontran)==1
        
        [~,colfinocup] = find (MTarryocupation==1);
        if  length (colfinocup)==1
            status = 0;
        end
        try
            if length(MTarryocupationtemp)>1
                nextposition =colfinocup(end) - colfinocup(end-1);
            else
                nextposition=1 ;
            end
        catch err
            return;
            %         nextposition=1
        end
        new_pos = pos(end) - (nextposition.*0.008);
        iarraysize = iarraysize-nextposition ;
        if nextposition>1
            testt = 1;
        end
        
        %     ocupationnumber(1)=sum(ocupationnumber(1:1+nextposition));
        %     ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
        ocupationnumber(1:positiontran+1)=ocupationnumber((1+nextposition):(positiontran+nextposition+1));
        %       ocupationnumber ((((positiontran-nextposition)+2):(positiontran+1))) = 0;
        MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
        
        %             if (length(ocupationnumber)-nextposition+1)>=1;
        %         ocupationnumber ((end-nextposition+1):end) = [];
        %             else
        %
        %                 ocupationnumber(1:end) = [];
        %                 break;
        %             end
        
        %   if lastnonzeromembranes(end)<lastnonzeroMT(end)
        %         testb =1
        %     end
    end
catch err
    if isempty(lastnonzeromembranes)==1
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber  = ocupationnumber ;
        MTarryocupation  = MTarryocupation ;
    end
    
end

%% Diffusion

%%Diffusion on membrane

if  transitionkind == 3 &&  positiontran~=0   && ocupationnumber (positiontran)~=0
    
    
    
    if  positiontran ~= lastnonzeromembranes(end) && ocupationnumber (positiontran+1)<=50
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        if positiontran~=1
            ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (3,:);
        else
            ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (3,:);
        end
        %     ocupationnumber (positiontran+1) =  ocupationnumber (positiontran+1)+1;
        %     ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
        %  if lastnonzeromembranes(end)<lastnonzeroMT(end)
        %         testb =1
        %     end
    end
    try
        if  positiontran == lastnonzeromembranes(end) && positiontran~=1&&lastnonzeromembranes(end)<=lastnonzeroMT(end) &&ocupationnumber (positiontran+1)<=50
            
            ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (3,:);
        end
    catch err
        return;
    end
    if positiontran == lastnonzeromembranes(end) && positiontran==1 &&lastnonzeromembranes(end)<=lastnonzeroMT(end) &&ocupationnumber (positiontran+1)<=50
        
        ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (3,:);
    end
    
    tranflag=3;
end

if  transitionkind == 4 &&  positiontran~=0   && ocupationnumber (positiontran)~=0
    
    if  positiontran~=1 && ocupationnumber (positiontran-1)<=50
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        %     ocupationnumber (positiontran-1) =  ocupationnumber (positiontran-1)+1;
        %     ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
        ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (4,:);
        %     if lastnonzeromembranes(end)<lastnonzeroMT(end)
        %         testb =1
        %     end
    end
    
    if  positiontran==1 && ocupationnumber (positiontran)<=50
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (4,:);
        %ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
    end
    tranflag=4;
end
%
%%Diffusion lattice microtubule membrane left rigth
if  transitionkind == 5 &&  positiontran~=0 && positiontran~=1 && MTarryocupation (positiontran)==1
    
    try
        if positiontran < lastnonzeromembranes(end)-1  && MTarryocupation (positiontran+1)==0
            new_pos = pos(end) + 0;
            iarraysize = iarraysize+0 ;
            MTarryocupation (positiontran+1) =  MTarryocupation (positiontran+1)+1;
            MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
        end
    catch err
        if isempty(lastnonzeromembranes)==1
            new_pos = pos(end) + 0;
            iarraysize = iarraysize+0 ;
            ocupationnumber  = ocupationnumber ;
            MTarryocupation  = MTarryocupation ;
        end
    end
    %
    
%     try
%     if positiontran == lastnonzeromembranes(end)-1  && MTarryocupation (positiontran+1)==0
%         new_pos = pos(end) + 0.008;
%       iarraysize = iarraysize+1 ;
%      ocupationnumber(2:positiontran+2)=ocupationnumber(1:positiontran+1);
%      MTarryocupation(2:positiontran+1)=MTarryocupation(1:positiontran);
%          MTarryocupation (positiontran+1) =  MTarryocupation (positiontran+1)+1;
%          MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
%     end
%     catch err
%       if isempty(lastnonzeromembranes)==1
%        new_pos = pos(end) + 0;
%      iarraysize = iarraysize+0 ;
%      ocupationnumber  = ocupationnumber ;
%     MTarryocupation  = MTarryocupation ;
%     end
%     end
    
    
end


if  transitionkind == 6 &&  positiontran~=0 && positiontran~=1 && MTarryocupation (positiontran)==1
    try
        if positiontran < lastnonzeromembranes(end)-1  && MTarryocupation (positiontran-1)==0
            new_pos = pos(end) + 0;
            iarraysize = iarraysize+0 ;
            MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
            MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
        end
    catch err
        if isempty(lastnonzeromembranes)==1
            new_pos = pos(end) + 0;
            iarraysize = iarraysize+0 ;
            ocupationnumber  = ocupationnumber ;
            MTarryocupation  = MTarryocupation ;
        end
    end
    %
%     try
%     if positiontran == lastnonzeromembranes(end)-1  && MTarryocupation (positiontran-1)==0
%            nextposition=1;
%           new_pos = pos(end) - (nextposition.*0.008);
%          iarraysize = iarraysize-nextposition ;
%     %
%     %
%          ocupationnumber(1)=sum(ocupationnumber(1:1+nextposition));
%          ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
%     ocupationnumber(1:positiontran+1)=ocupationnumber((1+nextposition):(positiontran+nextposition+1));
%     %
%            MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
%           MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
%     end
%     catch err
%       if isempty(lastnonzeromembranes)==1
%       new_pos = pos(end) + 0;
%     iarraysize = iarraysize+0 ;
%     ocupationnumber  = ocupationnumber ;
%     MTarryocupation  = MTarryocupation ;
%     end
%     end
    
    
%             if any (MTarryocupation(MTarryocupation~=0)>1)
%          testre = 1;
%     
%      end
%     if lastnonzeromembranes(end)<lastnonzeroMT(end)
%             testb =1
%         end
end



%   end

if MTarryocupation(1)==0
    MTarryocupation(1)=1;
    
end
if ocupationnumber(1)==0
    ocupationnumber(1)=1;
    
end
% test = find (MTarryocupation==-1);
%     if ~isempty(test)
%         return;
%     end

end

