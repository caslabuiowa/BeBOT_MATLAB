function [Cp, Transf] = HermitianInterpolation(x,tf)

%% x is a vector of initial/final conditions [x0 x0dot x0ddot ... xfddot xfdot xf]
% Cp is the Bernstein polynomials that satisfies these boundary conditions

if isrow(x)
    x = x';
end

n = length(x);

M = zeros(2,n);
M(1, n/2) = 1;
M(2, n/2 + 1) = 1;
Poly_dot = M*x;
A = tf/(length(Poly_dot)+2-1)*ones(length(Poly_dot),(length(Poly_dot)+2)/2)';
I = [tril(A,-1) ; -triu(A,1)];
E = DegElevMatrix(length(Poly_dot)-1,length(Poly_dot))';
midIdx = ceil((length(Poly_dot)+1)/2);
M_sel = eye(length(Poly_dot)+1); 
M_sel(midIdx, :) = []; 
    
B = zeros(length(Poly_dot) + 2, n);
B(1:(length(Poly_dot) + 2)/2, n/2 - 1) = 1;
B((length(Poly_dot) + 2)/2 + 1:end, n/2 + 1 + 1) = 1;

T = I*M_sel*E;

Transf = T*M + B;


for i = 2:length(x)/2 - 1

    A = tf/(2*i+1)*ones(2*i,(2*i+2)/2)';
    I = [tril(A,-1) ; -triu(A,i)];
    E = DegElevMatrix(2*i-1,2*i)';
    midIdx = ceil((2*i+1)/2);
    M_sel = eye(2*i+1); 
    M_sel(midIdx, :) = []; 
    
    B = zeros(2*i + 2, n);
    B(1:(2*i + 2)/2, n/2 - i) = 1;
    B((2*i + 2)/2 + 1:end, n/2 + 1 + i) = 1;

    T = I*M_sel*E;

    Transf = T*Transf + B;

end


Cp = Transf*x;


end


