function CM = DegReduction(CN,M)

N = length(CN) - 1;
r = N-M;

for k = 1 : M+1
    sum = 0
    for j = 1 : k
        sum = sum + (-1)^(k-j)*nchoosek_mod(k-j+r-1,r-1)*nchoosek_mod(N,j-1)/nchoosek_mod(N-r,k-1)*CN(j);
    end
    CM(k) = sum;
end


end









function binom = nchoosek_mod(N,k)
% This function produces the same output as MATLAB's nchoosek, but it's more efficient 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

binom = 1;
    for j = 1:k
        binom = binom*(N-(k-j));
        binom = binom/j;
    end
end