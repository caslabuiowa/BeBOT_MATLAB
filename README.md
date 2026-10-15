# BeBOT_MATLAB

This repository is the MATLAB version of BeBOT (https://github.com/caslabuiowa/BeBOT).

## Fast minimum distance

`BeBOT/MinDistBernstein2PolygonFast.m` and `BeBOT/MinDistBernstein2BernsteinFast.m` compute the minimum distance between a Bernstein polynomial and a convex polytope (or point), and between two Bernstein polynomials. They are drop-in alternatives to `MinDistBernstein2Polygon` and `MinDistBernstein2Bernstein`: a best-first branch and bound in which the GJK distance between convex hulls of the Bernstein coefficients gives a lower bound, points on the curve give an upper bound, and pieces are subdivided with the de Casteljau algorithm until the result is certified within the tolerance `epsilon` (default `1e-6`).

```matlab
[dist, t, pt, info] = MinDistBernstein2PolygonFast(cpts, poly, epsilon);
[dist, t1, t2, info] = MinDistBernstein2BernsteinFast(cpts1, cpts2, epsilon);
```

`cpts` are the Bernstein coefficients (one column per coefficient, time interval normalized to [0,1]) and `poly` the vertices of the polytope (one column per vertex). Helpers: `gjk_lite.m` (GJK distance with witness weights), `split_bern.m` (de Casteljau subdivision), `eval_bern.m` (evaluation). On the examples of Fig. 1 in the BeBOT book chapter, they are about 3–6x faster than the original functions in MATLAB R2024b.
