function x = eval_bern(C, s)
%EVAL_BERN Point of the Bernstein polynomial with coefficients C at s in [0,1]
%   (de Casteljau recursion, numerically stable).
N = size(C, 2) - 1;
for r = 1:N
    C(:,1:N+1-r) = (1-s)*C(:,1:N+1-r) + s*C(:,2:N+2-r);
end
x = C(:,1);
end
