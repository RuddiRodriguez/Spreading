% function []=get_statisctics_v2(orderpto,posinterp,posinterpmicro)

pks = 0 ;
locs = 0;
tube = 0 ;
small = 0;
nodeform=0;
vector =0 ;
% forcesq=0;
meanf=0;
meanf=0;
detach = 0;
lengtbedetach=NaN;
A=1:size(density,1);
B=[];
C=setdiff(A,B);
sliding=0;
trescri=0.3;
vector=0;
for j=C
 try   
    order_intert=orderpto{j, 1} ;
    order_intert=round (order_intert,1,'decimal');
    yyt = posinterp{j, 1} ;
    micropost = posinterpmicro{j,1};
    
    [pks,locs]= findpeaks (micropost,'MinPeakProminence',0.20);
    if isempty(pks);
        pks=micropost(end);
        locs = length(micropost);
    else
        pks = [pks';micropost(end)];
        locs = [locs';length(micropost)];
    end
    %    try
    %      figure ; plot (timeto(1:end),yyt,timeto(1:end),micropost,timeto(locs),micropost(locs),'o');hold on
    %    catch err
    %        figure ; plot (timeto(1:end-1),yyt,timeto(1:end-1),micropost,timeto(locs),micropost(locs),'o');hold on
    %    end
    %MinIdx=zeros(length(locs),1);
    % figure (3); plot (micropost);hold on ; plot (yyt); hold on ; plot(order_intert) ; hold off;
    % ylim([0 2])
    
    [pksindexrow,pksindexcol] = find((pks>trescri));
    
    locs = locs (pksindexrow);
    DataInv = 1.01*max(micropost ) - micropost ;
    [Minima,MinIdx] = findpeaks(DataInv,'MinPeakProminence',0.12);
    %     MinIdx(1,1)=1;
    %      MinIdx(end,1)=length(micropost);
    Minima = micropost (MinIdx);
    
    %     timespeaks= (timeto(locs));
    %     timespeaksmin= (timeto(MinIdx));
    if isempty(MinIdx);
        locs=locs(1);
    end
    %     try
    %         plot (timeto(1:end),yyt,timeto(1:end),micropost,timeto(locs),micropost(locs),'o',timeto(MinIdx),micropost(MinIdx),'s');hold off
    %     catch err
    %          plot (timeto(1:end-1),yyt,timeto(1:end-1),micropost,timeto(locs),micropost(locs),'o',timeto(MinIdx),micropost(MinIdx),'s');hold off
    %     end
    i=0;
     if length(locs)>=1
    for i=1:length(locs)
        if i ==1
            order_inter=order_intert(1:locs(1)) ;
            yy = yyt(1:locs(1)) ;
            micropos = micropost(1:locs(1));
                       %  figure(1);plot(yy);hold on ; plot (micropos);hold on;plot (order_inter,'o');hold off
            %             ylim([0 2])
            yy=yy-0.2;
            micropos =micropos-0.2;
        else
            
            order_inter=order_intert(MinIdx(i-1):locs(i)) ;
            yy = yyt(MinIdx(i-1):locs(i)) ;
            micropos = micropost(MinIdx(i-1):locs(i));
%                         figure;plot(yy);hold on ; plot (micropos);hold on;plot (order_inter,'o');hold off
%                         ylim([0 1])
%                         drawnow;
                        
                        
        end
        
        if length(order_inter)>5
            
            [lengtthmax]=(yy(find (order_inter>0.2)));
            [rowl,coll]=((find (order_inter>0.2)));
            [rowly,colly]=((find (yy>0.2)));
            [rowlm,collm]=((find (order_inter<0.2)));
            x=diff(collm);
            b=find([x inf]>1);
            c=diff([0 b]);
            d=cumsum(c);
            
            j;
            
            if ~isempty( collm)
                
                
                if ~isempty( b) && order_inter(end-2)>trescri && ~isempty( colly)
                    detach=[detach;1];
                    lengtbedetach = [lengtbedetach;yy(collm(d(end)))];
                    [fitresult, gof] = createFit_linear_sliding_check(1:length(yy(d:end)), yy(d:end));
                    if fitresult.p1>0.01 && fitresult.p1>0
                      sliding=[sliding;1];  
                    end
                else
                    detach=[detach;0];
                    lengtbedetach = [lengtbedetach;NaN];
                end
                
                if ~isempty( colly)
                    j;
                    %         yyle=yy(coll(1));
                    
                    
                    if ~isempty(coll)&& yy(collm(d(end)))<0.6
                        nodeform=[nodeform;1];
                        tube=[tube;0];
                        small=[small;0];
                    end
                    
                    if ~isempty(coll)&& yy(collm(d(end)))>0.6 %&&  pksm>1
                        tube=[tube;1];
                        small=[small;0];
                        nodeform=[nodeform;0];
%                         figure(1);plot(yy);hold on ; plot (micropos);hold on;plot (order_inter,'o');hold off
                    end
                    if ~isempty(coll) && yy(collm(d(end)))>0.2 &&  yy(collm(d(end)))<0.6
                        small=[small;1];
                        nodeform=[nodeform;0];
                        tube=[tube;0];
                        
                    end
                    
                    if isempty(coll)
                        tube=[tube;1];
                        small=[small;0];
                        nodeform=[nodeform;0];
                    end
                    
                else
                    tube=[tube;0];
                    small=[small;0];
                    nodeform=[nodeform;1];
                    
                    
                    
                end
                %    figure (1); plot (timeto,micropost,timeto,yyt)
                
            else
            tube=[tube;0];
            small=[small;0];
            nodeform=[nodeform;1];
            end
            
        else
            tube=[tube;NaN];
            small=[small;NaN];
            nodeform=[nodeform;NaN];
        end
        
    end
    
    
    
    
    
    
    %  if lengtthmax(1)<trescri
    %      nodeform=nodeform+1;
    %
    %  end
    j;
%     if i>1
    vector (j,1:8)=[length(tube)-1;sum(tube)+sum(small);sum(tube);sum(small);sum(nodeform);sum(detach);sum(sliding);i];
%     else
%      vector (j,1:7)=[NaN];   
%     end
    tube = 0;small=0;nodeform=0;timespeaks=0;detach=0;sliding=0;
 end
 catch
 end
end
% forcesq=6.28.*sqrt(vectorparameters(:,1).*2e-20)
% meanf=mean(forcesq).*1e12
% meanf=std(forcesq).*1e12

for i=1:size (density,1)
    data=vinterptoMT{i,1};
    data(data<=0)=NaN;
    vinterptoMTn{i,1}=data.*60;
    meanspped (i) = nanmean(data).*60;
    
    
    
end

[corow,~]=find (vector(:,1)==0);
vector(corow,:)=NaN;
meanspeedt = mean(meanspped);
stdspedtotal = (std (meanspped,[],2))./sqrt(length(meanspped));
meanprobtu=nanmean(vector(:,3)./vector (:,1));
meanprobde=nanmean(vector(:,6)./vector (:,1));
meanprobsliding=nanmean(vector(:,7)./vector (:,1));
lengtbedetach(lengtbedetach<=0)=NaN;
menalength = nanmean (lengtbedetach);

vectorre = [meanspeedt;meanprobtu;meanprobde;menalength;meanprobsliding;nansum(vector(:,8))];

vectorre=vectorre';