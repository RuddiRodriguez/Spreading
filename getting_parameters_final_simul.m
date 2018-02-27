
vectorparameters =zeros(length(parameterst),5);
for i=1:length(parameterst)
    
    [vectorparameters(i,:)] = parameterst{i,1}
%     sigma (i) = parameters{i,1}(2)
%     sigma (i) = parameters{i,1}(3)
%     sigma (i) = parameters{i,1}(4)
%     
end 
vectorparametersno= vectorparameters./(mean(vectorparameters));