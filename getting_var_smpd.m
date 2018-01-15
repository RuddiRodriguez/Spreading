function [post] = getting_var_smpd(pos)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pos1c = pos (1,1);
pos2c = pos (1,2);
pos3c = pos (1,3);

pos33c = [pos3c{1,1}];
pos22c = [pos2c{1,1}];
pos11c = [pos1c{1,1}];

if size (pos11c,2)~=1
    pos11c =pos11c';
end
if size (pos22c,2)~=1
    pos22c =pos22c';
end
if size (pos33c,2)~=1
    pos33c =pos33c';
end
try
post = [pos11c;pos22c;pos33c];
catch err
    post = [pos11c';pos22c';pos33c];
end
end

