function [new_pos,iarraysize,ocupationnumber,MTarryocupation,status] = unbound_transition(new_pos,positiontran,ocupationnumber,...
    pos,iarraysize,matrix_tmemla,matrix_tmemlap1,MTarryocupation,matrix_tMTla,matrix_tMTlap1,lastnonzeroMT,lastnonzeromembranes,status,MTarryocupationtemp)
%A bound transition from the membrane to the microtubule

%% unbound at the lattice

try
    if  positiontran~=0 && positiontran < lastnonzeroMT(end)  && MTarryocupation (positiontran)==1 && ocupationnumber (positiontran)<5
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        if positiontran~=1
            ocupationnumber (positiontran-1:positiontran+1) = ocupationnumber (positiontran-1:positiontran+1)+matrix_tmemla (1,:);
            MTarryocupation (positiontran-1:positiontran+1) = MTarryocupation (positiontran-1:positiontran+1)+matrix_tMTla (1,:);
        else
            ocupationnumber (positiontran:positiontran+1) = ocupationnumber (positiontran:positiontran+1)+matrix_tmemlap1 (1,:);
            MTarryocupation (positiontran:positiontran+1) = MTarryocupation (positiontran:positiontran+1)+matrix_tMTlap1 (1,:);
        end
    end
catch err
    if isempty(lastnonzeromembranes)==1
        new_pos = pos(end) + 0;
        iarraysize = iarraysize+0 ;
        ocupationnumber  = ocupationnumber ;
        MTarryocupation  = MTarryocupation ;
    end
    
end
%% unbound at the tip edge
try
    %unbound at the tip edge
    if positiontran == lastnonzeroMT(end) && MTarryocupation (positiontran)==1 && ocupationnumber (positiontran)<5
        
        [~,colfinocup] = find (MTarryocupation==1);
        if  length (colfinocup)==1
            status = 0;
             nextposition=colfinocup;
        end
        if length(MTarryocupation(MTarryocupation~=0))>1

            nextposition =colfinocup(end) - colfinocup(end-1);
        else
            nextposition=1 ;
        end
        new_pos = pos(end) - (nextposition.*0.008);
        iarraysize = iarraysize-nextposition ;
        
        MTarryocupation (positiontran) = MTarryocupation (positiontran)-1;
        ocupationnumber(positiontran)=ocupationnumber(positiontran)+1;
        ocupationnumber(1)=sum(ocupationnumber(1:1+nextposition));
        ocupationnumber(2:positiontran+1)=ocupationnumber((2+nextposition):(positiontran+nextposition+1));
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