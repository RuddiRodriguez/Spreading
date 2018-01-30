function [new_pos,iarraysize,ocupationnumber,MTarryocupation,controldensitynew] = diff_MT_right(new_pos,positiontran,ocupationnumber,...
    lastnonzeroMT,pos,iarraysize,matrix_tMTla,matrix_tMTlap1,lastnonzeromembranes,MTarryocupation,controldensitynew,MTarryocupationtemp)
%Diffusion to the right for the petides on the MT

%try
    
    if  positiontran~=0  && MTarryocupation (positiontran)==1 && MTarryocupation (positiontran+1)==0 && positiontran < lastnonzeroMT(end)
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        [controldensitynew] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
       %try
            if positiontran~=1
                MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (5,:);
                
            else
                MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (5,:);
                
            end
%             if positiontran==1
%             [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
%             
%             if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
%                 MTarryocupation =MTarryocupationtemp;
%                 
%             end
%             end
       % catch err
            if isempty(lastnonzeromembranes)==1
                new_pos = pos(end) + 0;
                iarraysize = iarraysize+0 ;
                ocupationnumber  = ocupationnumber ;
                MTarryocupation  = MTarryocupation ;
            end
       % end
        %
        
        %                 try
%                             if positiontran == lastnonzeroMT(end)  && MTarryocupation (positiontran+1)==0
%                                 new_pos = pos(end) + 0.008;
%                                 iarraysize = iarraysize+1 ;
%                                 ocupationnumber(2:positiontran+2)=ocupationnumber(1:positiontran+1);
%                                 MTarryocupation(2:positiontran+1)=MTarryocupation(1:positiontran);
%                                 MTarryocupation (positiontran+1) =  MTarryocupation (positiontran+1)+1;
%                                 MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
%                             end
        %                 catch err
        %                     if isempty(lastnonzeromembranes)==1
        %                         new_pos = pos(end) + 0;
        %                         iarraysize = iarraysize+0 ;
        %                         ocupationnumber  = ocupationnumber ;
        %                         MTarryocupation  = MTarryocupation ;
        %                     end
        %                 end
    end
% catch err
%     return;
% end
end