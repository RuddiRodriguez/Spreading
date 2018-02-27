
for j =1:size(ocupationnumbert)
    ocupationnumber=(ocupationnumbert{j, 1});
    MTarryocupation=(MTarryocupationt{j, 1});
for i =1:size(ocupationnumber,1)
    switch i
        case 1
    density (j,i) = (sum(ocupationnumber(i,2:end))+sum(MTarryocupation(i,2:end)))./0.240;
        otherwise
            density (j,i) = (sum(ocupationnumber(i,:))+sum(MTarryocupation(i,:)))./0.240;
        
    end
end
end