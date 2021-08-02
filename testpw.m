clear all
close all

% Cp = [0 1 2 0 -1 0 4 6 2 2 2 2];
% tknots = [0 2 5 6]'
% N = 3;
% t = linspace(0,6,1000);
% [tnodes,w,Dm] = CompBeBOT(N,tknots);
% CoBf = CompositeBernsteinPoly(Cp,tknots,t);
% plot(t,CoBf,'Color','g'); hold on
% plot(tnodes,Cp,'b*')



f = @(t) 1./(sqrt(1+(12-t).^2));
t = 0:0.01:12

N = 41;
tnodes = linspace(t(1),t(end),N+1);
Bf = BernsteinPoly(f(tnodes),t);

tknots = [0 8 11 12];
N = 5;
[tnodes,w,Dm] = PiecewiseBeBOT(N,tknots);
CoBf = PiecewiseBernsteinPoly(f(tnodes),tknots,t)


plot(t,f(t),'Color','b'); hold on
plot(t,Bf,'Color','r');
plot(t,CoBf,'Color','g');
plot(tnodes,f(tnodes),'b*')



% grid on

