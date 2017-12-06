A = zeros(1,4);
 for k = 1:7
   p = randperm(4); %  <-- This is the source of randomness
   
     A(1,p(ix)) = A(1,p(ix)) + 1;
   
 end
 A = [A(:,1),zeros(4,1),A(:,2:4)]; % Insert a column of four zeros
 A = [A(1:3,:);zeros(1,5);A(4,:)]; % Insert a row of five zeros
 A(1,1) = 0; A(3,3) = 0; A(5,5) = 0; % Set non-zero diagonal elements to 0