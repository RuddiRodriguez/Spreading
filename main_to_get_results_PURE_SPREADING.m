

try 
        F = uigetdir(dirtemp);
  catch err
      dirtemp = 'D:\Users_data\Ruddi\Dropbox\Data_analysis\Data_analysis\simulation\spreading\forces';
      F = uigetdir(dirtemp);
end
  
 dirtemp = fullfile (sprintf('%s',F));



  db = dir([dirtemp, '/**/*.mat']);
%db = dir([dirtemp, '/**/green/*.txt']);
A=1:size(db,1);
B=[];
C=setdiff(A,B);
% vectorfi=NaN(length(db),6);
% number=NaN(length(db),1);
% vectorforcet=NaN(length(db),3);
bin_means=NaN(length(db),12)
edges=NaN(length(db),12)
c=0;
for q =C;
    q
     path=fullfile(db(q).folder,db(q).name);
     
     
     load(sprintf('%s',path));
    
    Getting_results_from_spreading_to_plot_agains_displacement;
for cId = 1 : numel( vinterpto )
c=c+1
    vectorfiv(c,1:length(vinterpto{cId}))=vinterpto{cId}  ;
    vectorfipos(c,1:length(posttt{cId}))=posttt{cId}  ;
%     vectorbi = [posttt{cId};vinterpto{cId}]
%    [vectorress,edges(c,:),bin_means(c,:),bin_meansstd,N]= New_bining_data ( vectorbi',7,1,2)
 end
    
%      if vectorre(1,2)<=1
%     vectorfi(q,:)=vectorre;
%     number (q,:)=j;
%     vectorforce = [forcesq;meanprobtu;length(vectorparameters(:,1))];
%     vectorforcet(q,:)=vectorforce;
%      else
%          vectorfi(q,:)=NaN;
%      end
    clearvars -except db vectorfiv vectorfipos c;
end
 
% for l=1:64
%     l;
% vectorbi = [vectorfipos(l,:);vtrans(1,:)]
% [vectorress,edges1,bin_means1,bin_meansstd,N]= New_bining_data ( vectorbi',12,1,2);
% bin_means(l,1:length(bin_means1))=bin_means1;
% edges(l,1:length(edges1))=edges1;
% end

% [vectorressi,edges,bin_means,bin_meansstd,N]=New_bining_data ( vectorforcet,7,1,1);
% figure (1);errorbar (edges(1:end-1),bin_means,bin_meansstd./sqrt(N),'-o');hold on;
% [vectorresexp,edges,bin_means]=New_bining_data (TAC0v3,3,1,2);
% figure (1);plot (edges(1:end-1),bin_means,'-o');
% [vectorres,edges,bin_means,bin_meansstd,N]=New_bining_data (vectorfi0,10,1,2);
% figure (1);errorbar (edges(1:end-1),bin_means,bin_meansstd./sqrt(N),'-o');
% [vectorres,edges,bin_means,bin_meansstd,N]=New_bining_data (TAC200,3,2,5);
% figure (1);errorbar (edges(1:end-1),bin_means,bin_meansstd./sqrt(N),'-o');

