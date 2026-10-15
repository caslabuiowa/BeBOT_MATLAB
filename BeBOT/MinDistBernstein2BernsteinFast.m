function [dist, t1, t2, info] = MinDistBernstein2BernsteinFast(cpts1, cpts2, epsilon)
%MINDISTBERNSTEIN2BERNSTEINFAST Minimum distance between two Bernstein
%polynomials, within tolerance epsilon.
%   [dist, t1, t2] = MinDistBernstein2BernsteinFast(cpts1, cpts2, epsilon)
%   cpts1, cpts2: D x (N1+1), D x (N2+1) coefficients (D = 2 or 3), t in
%   [0,1]; epsilon: absolute tolerance (default 1e-6).
%   dist = ||x1(t1) - x2(t2)||, with dist - epsilon <= true minimum <= dist.
%
%   Best-first branch and bound (see MinDistBernstein2PolygonFast): lower
%   bounds are GJK distances between the convex hulls of the coefficients
%   of pairs of sub-curves; upper bounds are distances between actual
%   points of the curves (end points and GJK witness points); the pair
%   with the smallest lower bound is subdivided first, both curves at the
%   parameters of the GJK witness points.
if nargin < 3, epsilon = 1e-6; end
N1 = size(cpts1, 2) - 1;  N2 = size(cpts2, 2) - 1;  ngjk = 0;
% upper bound from the end points
E1 = cpts1(:,[1 end]);  E2 = cpts2(:,[1 end]);  alpha = inf;
for i = 1:2
    for j = 1:2
        dij = norm(E1(:,i) - E2(:,j));
        if dij < alpha, alpha = dij; t1 = i - 1; t2 = j - 1; end
    end
end
[lb, l1, l2] = gjk_lite(cpts1, cpts2);  ngjk = ngjk + 1;
C1 = {cpts1};  C2 = {cpts2};  ab = [0 1 0 1];  LB = lb;
S = [split_par(l1, N1), split_par(l2, N2)];
nsplit = 0;
while ~isempty(LB)
    [lbmin, k] = min(LB);
    if lbmin >= alpha - epsilon, break, end
    P = C1{k};  Q = C2{k};  abk = ab(k,:);  s = S(k,:);
    C1(k) = [];  C2(k) = [];  ab(k,:) = [];  LB(k) = [];  S(k,:) = [];
    [PA, PB, xp] = split_bern(P, s(1));  [QA, QB, xq] = split_bern(Q, s(2));
    nsplit = nsplit + 1;
    m1 = abk(1) + s(1)*(abk(2) - abk(1));  m2 = abk(3) + s(2)*(abk(4) - abk(3));
    % upper bound: split points and end points of the new pieces
    cand = [xp, xp, xp, P(:,1), P(:,end); xq, Q(:,1), Q(:,end), xq, xq];
    par  = [m1 m1 m1 abk(1) abk(2); m2 abk(3) abk(4) m2 m2];
    D = sqrt(sum((cand(1:end/2,:) - cand(end/2+1:end,:)).^2, 1));
    [dmin, c] = min(D);
    if dmin < alpha, alpha = dmin; t1 = par(1,c); t2 = par(2,c); end
    kids1 = {PA, [abk(1) m1]; PB, [m1 abk(2)]};
    kids2 = {QA, [abk(3) m2]; QB, [m2 abk(4)]};
    for a = 1:2
        for b = 1:2
            [lbc, la, lb2] = gjk_lite(kids1{a,1}, kids2{b,1});  ngjk = ngjk + 1;
            if lbc < alpha - epsilon
                % witness points of the sub-curves give an upper bound
                sa = (0:N1)*la/N1;  sb = (0:N2)*lb2/N2;
                xa = eval_bern(kids1{a,1}, sa);  xb = eval_bern(kids2{b,1}, sb);
                dw = norm(xa - xb);
                ia = kids1{a,2};  ib = kids2{b,2};
                if dw < alpha
                    alpha = dw;  t1 = ia(1) + sa*(ia(2) - ia(1));  t2 = ib(1) + sb*(ib(2) - ib(1));
                end
                if lbc < alpha - epsilon
                    C1{end+1} = kids1{a,1};  C2{end+1} = kids2{b,1}; %#ok<AGROW>
                    ab(end+1,:) = [ia ib];  LB(end+1) = lbc; %#ok<AGROW>
                    S(end+1,:) = [split_par(la, N1), split_par(lb2, N2)]; %#ok<AGROW>
                end
            end
        end
    end
end
dist = alpha;
info = struct('subdivisions', nsplit, 'gjk', ngjk);
end

function s = split_par(lam, N)
s = min(max((0:N)*lam/N, 0.1), 0.9);
end
