
data=0;
for t = 1:size(results,2)
   
    for tt=1:size(results,1)
%        data(tt,1:length(([results{tt,t,2}(:,3)])))=([results{tt,t,2}(:,3)]);
        plot([results{tt,t,1}(:,3)]);hold on;
    end
end