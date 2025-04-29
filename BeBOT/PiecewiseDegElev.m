function Cp_elev = PiecewiseDegElev(Cp,K,N,M)
% Cp = control points of piecewise Bernstein polynomial consisting of K polynomials of order N
% M is the desired order

Cp_m = reshape(Cp,N+1,K);

for i=1:K
    E = DegElevMatrix(N,M);
    Cp_elev(:,i) = Cp_m(:,i)'*E;
end

%Cp_elev = reshape(Cp_elev,K*(2*N+1),1);
Cp_elev = reshape(Cp_elev,[],1);
end









