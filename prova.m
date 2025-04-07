clear all
close all

Cp = [0 1 -1 2 2 1 0 1 1 3 3 3 3 1 0 -1];
N = 3;
tknots = [0 1 3 3.5 6];
t = linspace(0,tknots(end),1000);
lambda = 0.6;

[tnodes,~,~] = PiecewiseBeBOT(N,tknots);

figure 
plot(tnodes,Cp,'o','color','r') 
hold on
plot(t,PiecewiseBernsteinPoly(Cp,tknots,t),'color','r') 

[Cpout, Pos, tknotsout] = PiecewisedeCasteljau(Cp, tknots, lambda);

[tnodes,~,~] = PiecewiseBeBOT(N,tknotsout);


plot(tnodes,Cpout,'o','color','b') 
plot(t,PiecewiseBernsteinPoly(Cpout,tknotsout,t),'color','b') 
