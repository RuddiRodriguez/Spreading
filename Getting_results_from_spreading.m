 results = resultst
data=0;
pos=0;
times=0;
controldensity=0;
vinterpto = 0;
vinterptoM=0;
orderpto=0;
vinterpto = cell (size(results,1),1); 
 
for t = 1:size(results,1)
   
    for tt=1:size(results,3)
        pos=([results{t,1}(:,2)]);
        
        times=([results{t,1}(:,1)]);
          controldensity=([results{t,1}(:,3)]);
         
          
%          plot_mt_pos ( MTL, pos,times)
         plot(times,pos);hold on;
       [vector,yy,vinterp,timeto] = sampling_data(times,pos,controldensity,1,90);
      
        
          %vinterpto(t,1:length(vinterp))=vinterp;
            vinterpto{t}=vinterp;
 
        
    end
    
end
saving_to_txtfile;
function plot_mt_pos ( MTL, pos,times)
MTp=MTL(1:length(times)); 
memp = pos(1:length(times));
figure ; plot(times,MTp-memp);
figure
plot (times,MTL(1:length(times)),times,pos(1:length(times)))
    end