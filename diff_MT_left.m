function [new_pos,iarraysize,ocupationnumber,MTarryocupation] = diff_MT_left(new_pos,positiontran,ocupationnumber,...
    lastnonzeroMT,pos,iarraysize,lastnonzeromembranes,MTarryocupation,controldensitynew,MTarryocupationtemp)


if  positiontran~=0 && positiontran~=1 && MTarryocupation (positiontran)==1
    try
        if positiontran < lastnonzeroMT(end)
            new_pos = pos(end) + 0;
            iarraysize = iarraysize+0 ;
            if   positiontran~=2 && MTarryocupation (positiontran-1)==0
                MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
                MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
            end
            if  positiontran==2 && MTarryocupation (positiontran-1)<=50
                MTarryocupation (positiontran) = MTarryocupation (positiontran)-0;
                MTarryocupation (positiontran-1) =  MTarryocupation (positiontran-1)+1;
            end
            if  positiontran==1 
                MTarryocupation (positiontran) = MTarryocupation (positiontran);
                
            end
            
        end
%         [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
%         
%         if positiontran==2
%             if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
%                 MTarryocupation =MTarryocupation;
%             else
%                 MTarryocupation =MTarryocupation;
%             end
%             
%         end
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
    
    
end
end