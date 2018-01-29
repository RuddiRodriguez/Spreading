for i =1:size(ocupationnumber,1)
    switch i
        case 1
    density (i) = (sum(ocupationnumber(i,2:end))+sum(MTarryocupation(i,2:end)))./0.4;
        otherwise
            density (i) = (sum(ocupationnumber(i,:))+sum(MTarryocupation(i,:)))./0.4;
        
    end
end