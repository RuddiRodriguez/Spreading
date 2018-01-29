function [new_pos,iarraysize,ocupationnumber,tranflag] = diff_mem_right(new_pos,positiontran,ocupationnumber,...
    lastnonzeroMT,pos,iarraysize,matrix_tmemla,matrix_tmemlap1,lastnonzeromembranes,MTarryocupation,...
    controldensitynew,ocupationnumbertemp,tranflag)
%Diffusion to the right for the petides on the membrane


if  positiontran~=0   && ocupationnumber (positiontran)~=0
    if  positiontran ~= lastnonzeroMT(end)+1 && ocupationnumber (positiontran+1)<5
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        
        if positiontran~=1
            ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (3,:);
        else
            ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (3,:);
        end
    end
    if positiontran == lastnonzeroMT(end)+1  && lastnonzeromembranes(end)<=lastnonzeroMT(end) &&ocupationnumber (positiontran+1)<=50
        ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (3,:);
    end
%     if positiontran~=1
%     [controldensityneww] =density_calculation(MTarryocupation(2:end),ocupationnumber(2:end),new_pos);
%     if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
%         ocupationnumber =ocupationnumbertemp;
%     end
%     end
    tranflag=3;
end
end

