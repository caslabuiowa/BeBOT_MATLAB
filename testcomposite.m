close all

Cp = [0 -2 3 3 1 0 11 0 3 15];
tnodes = [0 3 5 7 10];
time = linspace(0,10,1000);

poly_t = CompositeBernsteinPoly(Cp,tnodes,time);

plot(linspace(tnodes(1),tnodes(end),length(poly_t)),poly_t);
grid on

