function Cp = PiecewiseBernsteinProduct(A,B,K,N)
% A is a piecewise Bernstein polynomial composed of K polynomials, each of
% order N
% B is a piecewise Bernstein polynomial composed of K polynomials, each of
% order N

A_m = reshape(A,N+1,K); % A_b is a matrix consisting of each individual Bernstein polynomial
B_m = reshape(B,N+1,K);
Cp_m = zeros(2*N+1,K);

for i = 1:K
    Cp_m(:,i) = BernsteinProduct(A_m(:,i)',B_m(:,i)');
end

Cp = reshape(Cp_m,K*(2*N+1),1);

end