function [new_pos,iarraysize,ocupationnumber,MTarryocupation] = bound_transition(new_pos,positiontran,ocupationnumber,...
    pos,iarraysize,matrix_tmemla,matrix_tmemlap1,MTarryocupation,matrix_tMTla,matrix_tMTlap1,lastnonzeroMT)
%A bound transition from the membrane to the microtubule

%% bound at the lattice


if positiontran~=0 && positiontran < lastnonzeroMT(end)+1 && ocupationnumber(positiontran)~=0 && MTarryocupation (positiontran)~=1
    new_pos = pos(end) + 0;
    iarraysize = iarraysize+0 ;
    if positiontran~=1
        ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (2,:);
        MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (2,:);
    else
        ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (2,:);
        MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (2,:);
    end
end
%% bound at the tip edge


if positiontran~=0 && positiontran == lastnonzeroMT(end)+1
    new_pos = pos(end) + 0.008;
    iarraysize = iarraysize+1 ;
    ocupationnumber(positiontran)=ocupationnumber(positiontran)-1;
    ocupationnumber(2:positiontran+1)=ocupationnumber(1:positiontran);
%     ocupationnumber(2)=randi([1 2],1,1);
    MTarryocupation (positiontran) = MTarryocupation (positiontran)+1;
    
end
end