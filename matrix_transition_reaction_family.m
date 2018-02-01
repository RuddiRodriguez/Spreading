function [matrix_tmemla,matrix_tMTla,matrix_tmemlap1,matrix_tMTlap1] = matrix_transition_reaction_family()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%In matrix_tmemlap1 the diffusion to the left in position 1 has reflecting
%boundary conditions
matrix_tmemla = [0 +1 0;0 -1 0 ; 0 -1 +1;+1 -1 0;0 0 0;0 0 0];
matrix_tmemlap1 = [ +1 0; -1 0 ;  -1 +1; +1 0; 0 0; 0 0];
matrix_tMTla = [0 -1 0;0 +1 0 ; 0 0 0;0 0 0;0 -1 +1;+1 -1 0];
matrix_tMTlap1 = [ -1 0; +1 0 ;  0 0; 0 0; -1 0; -1 0];
end

