function [dist, t, pt, info, trace] = MinDistBernstein2PolygonFast(cpts, poly, epsilon)
%MINDISTBERNSTEIN2POLYGONFAST Minimum distance between a Bernstein
%polynomial and a convex polytope (or a point), within tolerance epsilon.
%   [dist, t, pt] = MinDistBernstein2PolygonFast(cpts, poly, epsilon)
%   cpts: D x (N+1) Bernstein coefficients (D = 2 or 3), t in [0,1];
%   poly: D x M vertices of the polytope (M = 1 for a point);
%   epsilon: absolute tolerance (default 1e-6).
%   dist: distance of an actual point of the curve, x_N(t), from the
%   polytope, with dist - epsilon <= true minimum distance <= dist.
%   pt: closest point of the polytope.
%
%   Same principle as MinDistBernstein2Polygon (convex hull property, de
%   Casteljau subdivision, GJK), organized as a best-first branch and
%   bound: the sub-curve with the smallest lower bound (GJK distance
%   between its convex hull and the polytope) is subdivided first, at the
%   parameter of the GJK witness point; upper bounds are given by actual
%   points of the curve (end points and witness point of each sub-curve).
%   info: number of subdivisions and GJK calls.
%   trace: (optional) state of the algorithm at each iteration, used to
%   draw Figure 2 of the chapter.
if nargin < 3, epsilon = 1e-6; end
N = size(cpts, 2) - 1;  isPoint = size(poly, 2) == 1;
ngjk = 0;
    function [d, q] = dist_pt(x)        % distance of a curve point to the polytope
        if isPoint, q = poly; d = norm(x - poly);
        else, [d, ~, mu] = gjk_lite(x, poly); q = poly*mu; ngjk = ngjk + 1; end
    end
% upper bound from the two end points
[alpha, pt] = dist_pt(cpts(:,1));  t = 0;
[d1, q1] = dist_pt(cpts(:,end));
if d1 < alpha, alpha = d1; pt = q1; t = 1; end
% node list: coefficients, interval [a,b], lower bound, split parameter
[lb, lam] = gjk_lite(cpts, poly);  ngjk = ngjk + 1;
C = {cpts};  ab = [0 1];  LB = lb;  S = split_par(lam, N);
nsplit = 0;
doTrace = nargout > 4;  trace = {};  pruned = {};
if doTrace, trace{1} = snapshot(C, LB, ab, pruned, alpha, t, pt, [], cpts, poly); end
while ~isempty(LB)
    [lbmin, k] = min(LB);
    if lbmin >= alpha - epsilon, break, end
    Ck = C{k};  abk = ab(k,:);  s = S(k);
    C(k) = [];  ab(k,:) = [];  LB(k) = [];  S(k) = [];
    [CA, CB, x] = split_bern(Ck, s);  nsplit = nsplit + 1;
    tm = abk(1) + s*(abk(2) - abk(1));
    [dx, qx] = dist_pt(x);              % the split point is a point of the curve
    if dx < alpha, alpha = dx; pt = qx; t = tm; end
    kids = {CA, [abk(1) tm]; CB, [tm abk(2)]};
    for c = 1:2
        [lbc, lamc] = gjk_lite(kids{c,1}, poly);  ngjk = ngjk + 1;
        if lbc < alpha - epsilon
            C{end+1} = kids{c,1};  ab(end+1,:) = kids{c,2}; %#ok<AGROW>
            LB(end+1) = lbc;  S(end+1) = split_par(lamc, N); %#ok<AGROW>
        elseif doTrace
            pruned{end+1} = kids{c,1}; %#ok<AGROW>
        end
    end
    if doTrace
        % nodes whose lower bound is no longer below the incumbent are discarded
        drop = LB >= alpha - epsilon;
        pruned = [pruned, C(drop)];
        trace{end+1} = snapshot(C(~drop), LB(~drop), ab(~drop,:), pruned, alpha, t, pt, x, cpts, poly); %#ok<AGROW>
    end
end
dist = alpha;
info = struct('subdivisions', nsplit, 'gjk', ngjk);
end

function T = snapshot(C, LB, ab, pruned, alpha, t, pt, xsplit, cpts, poly)
% state of the branch and bound, for plotting
T.active = C;  T.lb = LB;  T.ab = ab;  T.pruned = pruned;
T.alpha = alpha;  T.t = t;  T.pt = pt;  T.xbest = eval_bern(cpts, t);
T.xsplit = xsplit;
T.witness = cell(size(C));            % closest points of the convex hulls
for k = 1:numel(C)
    [~, lam, mu] = gjk_lite(C{k}, poly);
    T.witness{k} = [C{k}*lam, poly*mu];
end
end

function s = split_par(lam, N)
% parameter of the GJK witness point (weights lam on the coefficients),
% kept away from the end points
s = min(max((0:N)*lam/N, 0.1), 0.9);
end
