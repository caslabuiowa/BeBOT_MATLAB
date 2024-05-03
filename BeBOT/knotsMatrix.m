function T = knotsMatrix(K,N)
 % Total number of elements in v
    totalElements = (N+1)*K;
    
    % Initialize T as a zero matrix
    T = zeros(K+1, totalElements);
    
    % Mark the first element of each segment
    for k = 1:K
        T(k, (k-1)*(N+1)+1) = 1;
    end
    
    T(end,end) = 1;

    T = T';
end

