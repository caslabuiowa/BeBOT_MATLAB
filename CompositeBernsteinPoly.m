function  poly_t  = CompositeBernsteinPoly(Cp,tnodes,time)

% INPUTS:
% 
% OUTPUT:
% 
% 
% Written by: Venanzio Cichella       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = length(tnodes)-1;
[dim, N] = size(Cp);
N = floor((N-1)/M);

poly_t = zeros(dim,length(time));

for i = 1 : M
    for k = 1 : length(time)   
        if time(k) >= tnodes(i) & time(k) <= tnodes(i+1) 
            t = time(k);
            if i < M
                poly_t(:,k) = BernsteinPoly(Cp((i-1)*N+1:i*N+1),t,tnodes(i),tnodes(i+1));
            else
                poly_t(:,k) = BernsteinPoly(Cp((i-1)*N+1:end),t,tnodes(i),tnodes(i+1));
            end
        end
    end
end


end







