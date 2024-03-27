function E = DegElevMatrixmod(N, M)
% INPUT: 
% N: Order of the curve to be elevated.
% M: Order to which it has to be elevated.
% OUTPUT: 
% Elevation matrix E.
% HOW TO USE: 
% Given a Bernstein poly of order N with control points cpN (row vector),
% the control points of the Bernstein poly degree elevated to order N,
% namely cpM (row vector), is given by cpM = cpN*E;

% Preallocate E with zeros for efficiency
E = zeros(M+1, N+1);

% Precompute all necessary binomial coefficients
binomN = arrayfun(@(k) nchoosek_mod(N, k), 0:N);
binomM = arrayfun(@(k) nchoosek_mod(M, k), 0:M);
binomR = arrayfun(@(k) nchoosek_mod(M-N, k), 0:(M-N));

% Vectorized filling of E
for i = 1:N+1
    for j = max(1, i):min(M+1, i+M-N)
        E(j,i) = binomN(i) * binomR(j-i+1) / binomM(j);
    end
end

% Transpose E to match the original function's output
E = E';
end








