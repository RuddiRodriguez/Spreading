function [post] = getting_var_smpd_8(pos)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
pos1c = pos (1,1);
pos2c = pos (1,2);
pos3c = pos (1,3);
pos4c = pos (1,4);
pos5c = pos (1,5);
pos6c = pos (1,6);
pos7c = pos (1,7);
pos8c = pos (1,8);
pos44c = [pos4c{1,1}];
pos33c = [pos3c{1,1}];
pos22c = [pos2c{1,1}];
pos11c = [pos1c{1,1}];
pos55c = [pos5c{1,1}];
pos66c = [pos6c{1,1}];
pos77c = [pos7c{1,1}];
pos88c = [pos8c{1,1}];

if size (pos11c,2)~=1
    pos11c =pos11c';
end
if size (pos22c,2)~=1
    pos22c =pos22c';
end
if size (pos33c,2)~=1
    pos33c =pos33c';
end
if size (pos44c,2)~=1
    pos44c =pos44c';
end
if size (pos55c,2)~=1
    pos55c =pos55c';
end
if size (pos66c,2)~=1
    pos66c =pos66c';
end
if size (pos77c,2)~=1
    pos77c =pos77c';
end
if size (pos88c,2)~=1
    pos88c =pos88c';
end

try
post = [pos11c;pos22c;pos33c;pos44c;pos55c;pos66c;pos77c;pos88c ];
catch err
    post = [pos11c';pos22c';pos33c];
end
end

