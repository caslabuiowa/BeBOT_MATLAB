function [d, lam, mu] = gjk_lite(P, Q)
%GJK_LITE Distance between the convex hulls of the columns of P and Q.
%   [d, lam, mu] = gjk_lite(P, Q) returns the Euclidean distance d between
%   conv(P) and conv(Q) (P, Q: D x m and D x n, D = 2 or 3), and convex
%   weights lam (m x 1) and mu (n x 1) of a pair of closest points,
%   P*lam and Q*mu. d = 0 if the hulls intersect.
%   Gilbert-Johnson-Keerthi algorithm on the Minkowski difference P - Q,
%   with the closest point of each simplex computed by enumerating its
%   faces (at most 15 for a tetrahedron). Lightweight replacement of
%   gjk.m (no input parsing, no rounding, no structs).
m = size(P, 2);  n = size(Q, 2);
I = 1;  J = 1;                          % simplex: pairs (I(k), J(k))
Wv = P(:,1) - Q(:,1);  a = 1;  v = Wv;
tol = 1e-12;
for it = 1:100
    vv = v'*v;
    if vv < tol^2, d = 0; break, end
    [~, i] = min(v'*P);  [~, j] = max(v'*Q);
    w = P(:,i) - Q(:,j);
    if vv - v'*w <= 1e-10*vv || any(I == i & J == j), break, end
    I(end+1) = i;  J(end+1) = j;  Wv(:,end+1) = w; %#ok<AGROW>
    [v, a] = closest_simplex(Wv);
    keep = a > 0;  I = I(keep);  J = J(keep);  Wv = Wv(:,keep);  a = a(keep);
    if numel(a) == size(P,1) + 1, v = zeros(size(v)); break, end   % origin inside
end
d = norm(v);
lam = zeros(m,1);  mu = zeros(n,1);
for k = 1:numel(a)
    lam(I(k)) = lam(I(k)) + a(k);  mu(J(k)) = mu(J(k)) + a(k);
end
end

function [v, a] = closest_simplex(Wv)
% closest point to the origin of conv(Wv), Wv: D x k, k <= 4
persistent SUB
if isempty(SUB)                 % non-empty subsets of {1,..,k}, k = 1..4
    SUB = cell(1,4);
    for kk = 1:4
        SUB{kk} = cell(1, 2^kk-1);
        for mask = 1:2^kk-1, SUB{kk}{mask} = find(bitget(mask, 1:kk)); end
    end
end
k = size(Wv, 2);  best = inf;  v = Wv(:,1);  a = [1; zeros(k-1,1)];
subs = SUB{k};
for q = 1:numel(subs)
    S = subs{q};
    w0 = Wv(:,S(1));
    if numel(S) == 1
        b = 1;  p = w0;
    else
        E = Wv(:,S(2:end)) - w0;  G = E'*E;
        if rcond(G) < 1e-14, continue, end
        c = -G\(E'*w0);  b = [1 - sum(c); c];
        if any(b < -1e-12), continue, end
        p = Wv(:,S)*b;
    end
    pp = p'*p;
    if pp < best - 1e-15
        best = pp;  v = p;  a = zeros(k,1);  a(S) = b;
    end
end
a = max(a, 0);  a = a/sum(a);
end
