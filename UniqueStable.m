function [AA, AI, BI] = UniqueStable(A)
% Equivalent to: unique(A(:), 'stable')
% Inptus: single, double (U)INT8/16/31/64, char
% Author: Jan Simon, Heidelberg 2017, License: CC BY-SA 3.0
nA = numel(A);
if nA == 0
   AA = [];
   AI = [];
   BI = [];
else
   [As, SV] = sort(A(:));
   if nargout == 1
      UV(SV) = ([true; diff(As) ~= 0]);
      AA     = A(UV);
     else  % Indices requested:
        UV      = ([true; diff(As) ~= 0]);
        UVs(SV) = UV;
        AI      = find(UVs).';
        AA      = A(UVs);
        % Complex creation of BI such that AA(BI) == A:
        v      = zeros(nA, 1);
        v(AI)  = 1:length(AI);    % Sequence related to AA
        vs     = v(SV);           % Sorted like A
        vf     = vs(vs ~= 0);     % Just the filled entries
        BI     = zeros(nA, 1);    % Pre-allocate
        BI(SV) = vf(cumsum(UV));  % Inflate multiple elements
     end
  end