function [arrayrates] = arrayrates_values_family_reaction(count,MTarryocupationtemp,numberpb,iMTLsize,MTarryocupation,ocupationnumber,rates,koof,arrayrates,positiontran,tranflag,densitylb,densitylu,v)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
ratestemp34 =(densitylu.*v);
ratestemp56 =(densitylb.*v);

if count==1 
        
        arrayrates (2,1:iMTLsize) = rates(1,2).*ocupationnumber;
        
        arrayrates (1,1:iMTLsize) = rates(1,1).*MTarryocupation;
         arrayrates (3,2:iMTLsize) = rates(1,3).*ocupationnumber(2:end);
        
% %         arrayrates (4,1:iMTLsize) = rates(1,4).*ocupationnumber;
% %         arrayrates (5,1:iMTLsize) = rates(1,5).*MTarryocupation;
           arrayrates (4,2:iMTLsize) = rates(1,4).*ocupationnumber(2:end);
          
           if size (arrayrates,1)==6
           arrayrates (5,2:iMTLsize) = rates(1,5).*MTarryocupation(2:end);
           arrayrates (6,2:iMTLsize) = rates(1,6).*MTarryocupation(2:end);
           arrayrates (5,1) = ratestemp56.*MTarryocupation(1);
           arrayrates (6,1) = ratestemp56.*MTarryocupation(1);
           end
 end
    
    if length (MTarryocupationtemp)< numberpb || length (MTarryocupationtemp) == 4
        
        arrayrates (2,1:iMTLsize) = rates(1,2).*ocupationnumber;
        
        arrayrates (1,1:iMTLsize) = rates(1,1).*MTarryocupation;
         arrayrates (3,2:iMTLsize) = rates(1,3).*ocupationnumber(2:end);
          
% %         arrayrates (4,1:iMTLsize) = rates(1,4).*ocupationnumber;
% %         arrayrates (5,1:iMTLsize) = rates(1,5).*MTarryocupation;
           arrayrates (4,2:iMTLsize) = rates(1,4).*ocupationnumber(2:end);
           
           if size (arrayrates,1)==6
           arrayrates (5,2:iMTLsize) = rates(1,5).*MTarryocupation(2:end);
           arrayrates (6,2:iMTLsize) = rates(1,6).*MTarryocupation(2:end);
           rrayrates (5,1) = ratestemp56.*MTarryocupation(1);
           arrayrates (6,1) = ratestemp56.*MTarryocupation(1);
           end
    end
    
    
    
    
    try
    if length (MTarryocupationtemp)> numberpb && positiontran~=1
        
        
        if tranflag==3 || tranflag==3
             arrayrates (2,positiontran-1:positiontran+1) = rates(1,2).*ocupationnumber(positiontran-1:positiontran+1);
              arrayrates (3,positiontran-1:positiontran+1) = rates(1,3).*ocupationnumber(positiontran-1:positiontran+1);
              arrayrates (4,positiontran-1:positiontran+1) = rates(1,4).*ocupationnumber(positiontran-1:positiontran+1);
        else
        
        [~,colocuone] = (find (MTarryocupation==1));
%         colocuone = myFind( MTarryocupation );
        %         try
        arrayrates (2,positiontran-1:positiontran+1) = rates(1,2).*ocupationnumber(positiontran-1:positiontran+1);
        %         catch err
        %             arrayrates=1;
        %         end
        %arrayrates (2,((colocuone(end)-numberpb+1):colocuone(end))) = rates(1,2).*MTarryocupation(((colocuone(end)-numberpb+1):colocuone(end)));
        
        arrayrates (1,1:(colocuone(end)-numberpb)) = koof.*MTarryocupation(1:(colocuone(end)-numberpb));
        arrayrates (1,1:iMTLsize) = rates(1,1).*MTarryocupation;
        arrayrates (3,positiontran-1:positiontran+1) = rates(1,3).*ocupationnumber(positiontran-1:positiontran+1);
        % arrayrates (4,1:iMTLsize) = rates(1,4).*ocupationnumber;
%         arrayrates (5,1:iMTLsize) = rates(1,5).*MTarryocupation;
         arrayrates (4,positiontran-1:positiontran+1) = rates(1,4).*ocupationnumber(positiontran-1:positiontran+1);
         if size (arrayrates,1)==6
         arrayrates (5,positiontran-1:positiontran+1) = rates(1,5).*MTarryocupation(positiontran-1:positiontran+1);
         arrayrates (6,positiontran-1:positiontran+1) = rates(1,6).*MTarryocupation(positiontran-1:positiontran+1);
         end
        end
    end
    catch err
        return;
    end
    if length (MTarryocupationtemp)> numberpb && positiontran==1
        
        [~,colocuone] = (find (MTarryocupation==1));
%         colocuone = myFind( MTarryocupation );
        %         try
        arrayrates (2,1:iMTLsize) = rates(1,2).*ocupationnumber;
        %         catch err
        %             arrayrates=1;
        %         end
        %arrayrates (2,((colocuone(end)-numberpb+1):colocuone(end))) = rates(1,2).*MTarryocupation(((colocuone(end)-numberpb+1):colocuone(end)));
        try
        arrayrates (1,1:(colocuone(end)-numberpb)) = koof.*MTarryocupation(1:(colocuone(end)-numberpb));
        arrayrates (1,1:iMTLsize) = rates(1,1).*MTarryocupation;
        arrayrates (3,positiontran:positiontran+1) = rates(1,3).*ocupationnumber(positiontran:positiontran+1);
        catch err
            return;
        end
        % arrayrates (4,1:iMTLsize) = rates(1,4).*ocupationnumber;
%         arrayrates (5,1:iMTLsize) = rates(1,5).*MTarryocupation;
         arrayrates (4,positiontran:positiontran+1) = rates(1,4).*ocupationnumber(positiontran:positiontran+1);
  
         if size (arrayrates,1)==6
         arrayrates (5,positiontran:positiontran+1) = rates(1,5).*MTarryocupation(positiontran:positiontran+1);
         arrayrates (6,positiontran:positiontran+1) = rates(1,6).*MTarryocupation(positiontran:positiontran+1);
         
         arrayrates (5,positiontran) = ratestemp56.*MTarryocupation(positiontran);
         arrayrates (6,positiontran) = ratestemp56.*MTarryocupation(positiontran);
         end
    end
arrayrates (3,1) = ratestemp34.*ocupationnumber(1);
arrayrates (4,1) = ratestemp34.*ocupationnumber(1);
if size (arrayrates,1)==6
    arrayrates (5,1) = ratestemp56.*MTarryocupation(1);
         arrayrates (6,1) = ratestemp56.*MTarryocupation(1);
       
end
end
