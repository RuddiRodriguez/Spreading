function [numberbtemp] = random_fix_sum(numberb,controlsize,density )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for hh =1:500
    numberbtemp = 0;
    sumtemp=0;
    h=0;
while (sumtemp)<numberb || length(numberbtemp)<numberb
    h=h+1;
    numberbtemp (h) = randi ([0 density],1,1); 
    sumtemp = sum  (numberbtemp);
    if sumtemp>numberb
        numberbtemp = numberbtemp(1:end-1);
    end
    if length (numberbtemp)>=controlsize   
     break;
 end 
end                                        
 
end
end

