function [new_pos,iarraysize,ocupationnumber,tranflag] = diff_mem_left(new_pos,positiontran,ocupationnumber,...
    pos,iarraysize,matrix_tmemla,matrix_tmemlap1,MTarryocupation,controldensitynew,ocupationnumbertemp,tranflag)
%Diffusion to the left of the petides on the membrane
% If the position is 1 everything is the same
% If the position is 2 we have reflecting boundary conditions the position
% 1 increase in 1 and the occupation number at the position 2 remain
% equally .However we still have a control parameter  ocupationnumbertemp that keep
%the density constant if we want

if  positiontran~=0   && ocupationnumber (positiontran)~=0
    if  positiontran~=1 && ocupationnumber (positiontran-1)<5 && positiontran~=2
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        try
        ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (4,:);
         catch err
            return
        end
        if positiontran==2 && ocupationnumber (positiontran-1)<5
            new_pos = pos(end) + 0;
            iarraysize = iarraysize+0 ;
            ocupationnumber(positiontran-1:positiontran) = ocupationnumber (positiontran-1:positiontran)+matrix_tmemlap1 (4,:);
        end
   
    if  positiontran==1 
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber = ocupationnumbertemp;%ocupationnumber (positiontran-1:positiontran)+matrix_tmemlap1 (4,:);
    end
%     if positiontran==2
%         [controldensityneww] =density_calculation(MTarryocupation(2:end),ocupationnumber(2:end),new_pos);
%         if controldensityneww <=(controldensitynew-5) || controldensityneww >=(controldensitynew+5)
%             ocupationnumber =ocupationnumbertemp;%here it is possible to use a constrain with ocupationnumbertemp to control the density
%         
%         end
        
        
    %end
    tranflag=4;
end
end
