function [A, B, x] = split_bern(C, s)
%SPLIT_BERN de Casteljau subdivision of a Bernstein polynomial at s in (0,1).
%   C: D x (N+1) coefficients; A, B: coefficients of the two pieces on
%   [0,s] and [s,1] (rescaled to [0,1]); x: point of the polynomial at s.
N = size(C, 2) - 1;  A = C;  B = C;  T = C;
for r = 1:N
    T(:,1:N+1-r) = (1-s)*T(:,1:N+1-r) + s*T(:,2:N+2-r);
    A(:,r+1) = T(:,1);  B(:,N+1-r) = T(:,N+1-r);
end
x = A(:,end);
end
