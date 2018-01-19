count =0;
xnew =0:0.1:1;
positiont = 0;
position=position1;
for i= 2:size(position,2)
    count = count +1
    [row, ~] = find(~isnan(position(:,i)));
    timen=position (1:length(row),1);
    timenn = timen./timen(end);
    positionnew = position (1:length(row),i);
    
    interponew = interp1(timenn,positionnew,xnew);
     figure(1) ; plot (xnew, interponew);hold on;
    
    positiont(1:length(xnew),count)=interponew;
    
end 
menaposition = nanmean(positiont,2);
stdposition = nanstd(positiont,[],2);
figure(2) ;plot (xnew,menaposition)