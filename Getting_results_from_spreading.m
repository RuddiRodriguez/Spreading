
data=0;
for t = 1:size(results,1)
   
    for tt=1:size(results,3)
        %data(tt,1:length(([results{tt,t,2}(:,3)])))=([results{tt,t,2}(:,3)]);
        plot([results{t,1}(:,end)]);hold on;
    end
end