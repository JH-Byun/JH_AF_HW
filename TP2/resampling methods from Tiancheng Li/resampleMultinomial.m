function [ indx ] = resampleMultinomial(w, N)
% [ indx ] = resampleMultinomial(w, N)
% Multinomial resampling method for particle filtering. 
% Author: Tiancheng Li,Ref:
% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% submit to IEEE Signal Processing Magazine, August 2013

% Input:
%       w    the input weight sequence 
%       N    the desired length of the output sequence(i.e. the desired number of resampled particles)
% Output:
%       indx the resampled index according to the weight sequence

if nargin == 1
  N = length(w);
end
w = w / sum(w);
Q = cumsum(w);

indx = zeros(1, N);
i = 0;
while i < N
    i = i + 1;
    sampl = rand;  % ~(0,1]
    j = 1;
    while Q(j) < sampl
        j = j + 1;
    end;
    indx(i)= j;
end
