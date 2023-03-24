function Prod = ProdMatrix(N)
%This function produces a matrix which can be used to compute ||x dot x||^2
% i.e. xaug = x'*x;
% xaug = reshape(xaug',[length(x)^2,1]);
% y = Square*xaug;
% or simply norm_square(x)
% This gives the same results as but is more efficient than BernsteinProduct(x,x)


T = zeros(2*N+1,(N+1)^2);

for j = 0:2*N
for i = max(0,j-N): min(N,j)
   if N >= i && N >= j-i && 2*N >= j && j-i >= 0
       T(j+1,N*i+j+1) = nchoosek_mod(N,i)*nchoosek_mod(N,j-i)/nchoosek_mod(2*N,j);
   end
end   
end

Prod = T;


end

