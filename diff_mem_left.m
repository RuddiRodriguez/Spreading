function [new_pos,iarraysize,ocupationnumber,tranflag] = diff_mem_left(new_pos,positiontran,ocupationnumber,...
    pos,iarraysize,matrix_tmemla,matrix_tmemlap1,MTarryocupation,controldensitynew,ocupationnumbertemp,tranflag)
%Diffusion to the left for the petides on the membrane

if  positiontran~=0   && ocupationnumber (positiontran)~=0
    if  positiontran~=1 && ocupationnumber (positiontran-1)<5 && positiontran~=2
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (4,:);
        if positiontran==2
        end
    end
    if  positiontran==1 && ocupationnumber (positiontran)<=5
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (4,:);
    end
    [controldensityneww] =density_calculation(MTarryocupation,ocupationnumber,new_pos);
    if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
        ocupationnumber =ocupationnumbertemp;
    end
    tranflag=4;
end
end
