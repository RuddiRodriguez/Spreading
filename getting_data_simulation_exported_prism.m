

index= find(isnan(density4(:,1))==1);
indext=0;
posititotal=0;
 figure;
for i= 1:length(index)
    
    time = density4(indext+1:index(i)-1,1);
    positi = density4(indext+1:index(i)-1,2);
    density = density4(indext+1:index(i)-1,3);
    timeinter = 0:0.1:max(time);
    yypositi = (interp1(time,positi(1:length(time)),timeinter))-0.2;
    plot (time,positi);hold on;plot(timeinter,yypositi,'ro');hold on
    indext= index(i);
    posititotal (1:length(yypositi),i)=yypositi;
    
end
 posititotal( posititotal==0)=NaN;
  meanpositi= nanmean(posititotal,2);
  stdpositi= nanstd(posititotal,[],2);
  vector = [timeinter' meanpositi stdpositi];
  figure (1); errorbar(timeinter,meanpositi,stdpositi);