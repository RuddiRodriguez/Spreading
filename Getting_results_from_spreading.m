
data=0;
pos=0;
times=0;
controldensity=0;
vinterpto = 0;
for t = 1:size(results,1)
   
    for tt=1:size(results,3)
        pos=([results{t,1}(:,2)]);
        
        times=([results{t,1}(:,1)]); 
         controldensity=([results{t,1}(:,end)]); 
        plot(times,pos-pos(1));hold on;
         [vector,yy,vinterp] = sampling_data(times,pos-pos(1),controldensity,1,20)
         vinterpto(t,1:length(vinterp))=vinterp;
    end
end