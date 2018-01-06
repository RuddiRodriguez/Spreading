function [ ocupationnumber,MTarryocupation,new_pos,pos,iarraysize,status,tranflag,controldensitynew] = transitions_no_EB_infinite_family_reaction_backup_6reaction(transitionkind,positiontran,ocupationnumber,MTarryocupation,MTarryocupationtemp,new_pos,pos,iarraysize,...
    matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1,ouft,inft,controldensitynew)
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

postemp = pos(end);
iarraysizetemp = iarraysize;
ocupationnumbertemp  = ocupationnumber;
MTarryocupationtemp = MTarryocupation;
ouft=ouft+0; 

switch transitionkind
    
    %% Diffusion
    case 3
        [new_pos,iarraysize,ocupationnumber,tranflag] = diff_mem_right(new_pos,positiontran,ocupationnumber,...
    lastnonzeroMT,pos,iarraysize,matrix_tmemla,matrix_tmemlap1,lastnonzeromembranes,MTarryocupation,...
    controldensitynew,ocupationnumbertemp,tranflag);
        %%Diffusion on membrane
        % Diffusion to the right
       
    case 4
        % Diffusion to the left
        [new_pos,iarraysize,ocupationnumber,tranflag] = diff_mem_left(new_pos,positiontran,ocupationnumber,...
    pos,iarraysize,matrix_tmemla,matrix_tmemlap1,MTarryocupation,controldensitynew,ocupationnumbertemp,tranflag);
        
        
        
    case 2
        [new_pos,iarraysize,ocupationnumber,MTarryocupation] = bound_transition(new_pos,positiontran,ocupationnumber,...
    pos,iarraysize,matrix_tmemla,matrix_tmemlap1,MTarryocupation,matrix_tMTla,matrix_tMTlap1,lastnonzeroMT);
        
%         %bound at the lattice
%         if transitionkind == 2 &&     positiontran~=0 && positiontran < lastnonzeroMT(end)+1 && ocupationnumber(positiontran)~=0 && MTarryocupation (positiontran)~=1
%             new_pos = pos(end) + 0;
%             iarraysize = iarraysize+0 ;
%             %     ocupationnumber (positiontran) = 1;
%             %     ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
%             if positiontran~=1
%                 ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (2,:);
%                 MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (2,:);
%             else
%                 ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (2,:);
%                 MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (2,:);
%             end
%             %     if ocupationnumber (positiontran)==0
%             %         ocupationnumber (positiontran)=randi ([1 3],1,1);
%             %     end
%             % MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (2,:);
%             % MTarryocupation (positiontran) = MTarryocupation (positiontran)+1;
%             
%             %  if lastnonzeromembranes(end)<lastnonzeroMT(end)
%             %         testb =1
%             %     end
%             
%         end
%         %bound at the tip edge
%         if transitionkind == 2 &&   positiontran~=0 && positiontran == lastnonzeroMT(end)+1 %&& MTarryocupation (positiontran)==0
%             
%             new_pos = pos(end) + 0.008;
%             iarraysize = iarraysize+1 ;
%             ocupationnumber(positiontran)=ocupationnumber(positiontran)-1;
%             ocupationnumber(3:positiontran+1)=ocupationnumber(2:positiontran);
%             ocupationnumber(2)=1;
%             %     ocupationnumber (positiontran) = ocupationnumber (positiontran)-1;
%             %     ocupationnumber (positiontran+1) = randi ([1 1],1,1);
%             %     if ocupationnumber (positiontran)==0
%             %         ocupationnumber (positiontran)=randi ([1 3],1,1);
%             %     end
%             
%             MTarryocupation (positiontran) = MTarryocupation (positiontran)+1;
%             
%             
%             %   if lastnonzeromembranes(end)<lastnonzeroMT(end)
%             %         testb =1
%             %     end
%             
%         end
        
    case 1
        
        
        %unbound at the lattice
        try
            if  transitionkind == 1 &&    positiontran~=0 && positiontran < lastnonzeroMT(end)  && MTarryocupation (positiontran)==1 && ocupationnumber (positiontran)<5
                
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
        
        
        
        
        
        
        
        % if positiontran == length ((ocupationnumber(ocupationnumber~=0)))-1 &&transitionkind == 1 && MTarryocupation (positiontran-1)==1
        %     return ;
        % end
        try
            %unbound at the tip edge
            if transitionkind == 1 &&  positiontran == lastnonzeroMT(end) && MTarryocupation (positiontran)==1 && ocupationnumber (positiontran)<5
                controltrnasitionlattice=0;
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
                MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
                ocupationnumber(positiontran)=ocupationnumber(positiontran)+1;
                ocupationnumber(1)=sum(ocupationnumber(1:1+nextposition));
                %     ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
                ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
                %       ocupationnumber ((((positiontran-nextposition)+2):(positiontran+1))) = 0;
                
                
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
            controltrnasitionlattice=0;
        end
        
        
    case 5
        try
            %
            %Diffusion lattice microtubule membrane left rigth
            %Diffusion to the rigth on the microtubule
            if  transitionkind == 5 &&  positiontran~=0  && MTarryocupation (positiontran)==1 && MTarryocupation (positiontran+1)==0 && positiontran < lastnonzeroMT(end)
                [controldensitynew] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
                try
                    %         if positiontran < lastnonzeroMT(end)  && MTarryocupation (positiontran+1)==0
                    %             new_pos = pos(end) + 0;
                    %             iarraysize = iarraysize+0 ;
                    %             MTarryocupation (positiontran+1) =  MTarryocupation (positiontran+1)+1;
                    %             MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
                    %         end
                    
                    if positiontran~=1
                        MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (5,:);
                        
                    else
                        MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (5,:);
                        inft=inft+1;
                    end
                    
                    [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
                    
                    if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
                        MTarryocupation =MTarryocupationtemp;
                        
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
                
                try
                    if positiontran == lastnonzeroMT(end)  && MTarryocupation (positiontran+1)==0
                        new_pos = pos(end) + 0.008;
                        iarraysize = iarraysize+1 ;
                        ocupationnumber(2:positiontran+2)=ocupationnumber(1:positiontran+1);
                        MTarryocupation(2:positiontran+1)=MTarryocupation(1:positiontran);
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
                
                
            end
        catch err
            return;
        end
        
    otherwise
        %Diffussion to the left on the microtubule
        if  transitionkind == 6 &&  positiontran~=0 && positiontran~=1 && MTarryocupation (positiontran)==1
            try
                if positiontran < lastnonzeroMT(end)
                    new_pos = pos(end) + 0;
                    iarraysize = iarraysize+0 ;
                    if   positiontran~=2 && MTarryocupation (positiontran-1)==0
                        MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
                        MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
                    end
                    if  positiontran==2 && MTarryocupation (positiontran-1)<=5
                        MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
                        MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
                    end
                end
                
                [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
                
                if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
                    MTarryocupation =MTarryocupationtemp;
                    
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
            %     if positiontran == lastnonzeroMT(end) % && MTarryocupation (positiontran-1)==0
            %
            %         [~,colfinocup] = find (MTarryocupation==1);
            %         if  length (colfinocup)==1
            %             status = 0;
            %         end
            %         try
            %             if length(MTarryocupationtemp)>1
            %                 nextposition =colfinocup(end) - colfinocup(end-1);
            %             else
            %                 nextposition=1 ;
            %             end
            %         catch err
            %             return;
            %             %         nextposition=1
            %         end
            %
            % new_pos = pos(end) - (nextposition.*0.008);
            %         iarraysize = iarraysize-nextposition ;
            %
            %         MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
            % ocupationnumber(positiontran+1)=ocupationnumber(positiontran+1)+1;
            %              ocupationnumber(1)=sum(ocupationnumber(1:1+nextposition));
            %         %     ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
            %         ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
            %
            %
            %
            %
            %
            %
            % %
            % %
            % %           new_pos = pos(end) - (nextposition.*0.008);
            % %          iarraysize = iarraysize-nextposition ;
            % %     %
            % %     %
            % %          ocupationnumber(1)=sum(ocupationnumber(1:1+nextposition));
            % %          ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
            % %     ocupationnumber(1:positiontran+1)=ocupationnumber((1+nextposition):(positiontran+nextposition+1));
            % %     %
            % %            MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
            % %           MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
            %     end
            %     catch err
            %       if isempty(lastnonzeromembranes)==1
            %       new_pos = pos(end) + 0;
            %     iarraysize = iarraysize+0 ;
            %     ocupationnumber  = ocupationnumber ;
            %     MTarryocupation  = MTarryocupation ;
            %     end
            %     end
            %
            
            if any (MTarryocupation(MTarryocupation~=0)>1)
                testre = 1;
                
            end
            if lastnonzeromembranes(end)<lastnonzeroMT(end)
                testb =1;
            end
        end
        
end

%   end

if MTarryocupation(1)~=1
    MTarryocupation(1)=1;
    
end
if any (MTarryocupation<0)
    testtt=1;
    
end


if ocupationnumber(1)~=1
    ocupationnumber(1)=1;
    
end
% [controldensitynew] =density_calculation(MTarryocupation,ocupationnumber,pos(end));
% if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
%                 ocupationnumber =ocupationnumbertemp;
%                  MTarryocupation =MTarryocupationtemp;
%                  new_pos = postemp;
%                 iarraysize = iarraysizetemp ;
%             end
% if ocupationnumber(2)~=0
%     ocupationnumber(2)=1;
%
% end
% test = find (MTarryocupation==-1);
%     if ~isempty(test)
%         return;
%     end
% MTocutemp = MTarryocupation(2:end);
% memocutemp = ocupationnumber(2:end);
% controldensityneww=(sum(MTocutemp(MTocutemp~=0))+sum(memocutemp(memocutemp~=0)))/new_pos ;
% if controldensityneww <controldensitynew %&& transitionkind==3
%
%     ctr=ctr+1;
%     transitionkindtest (ctr) =  transitionkind;
%     positiontrantest (ctr) = positiontran;
%
%     testt = 1;
% else
%     ctr=ctr+1;
%     transitionkindtest (ctr) = NaN;
%     positiontrantest (ctr) = NaN;
% end
end

