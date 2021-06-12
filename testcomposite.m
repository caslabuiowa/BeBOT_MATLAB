clear all
close all



f = @(t) 1./(sqrt(1+(12-t).^2));
t = 0:0.01:12

N = 41;
tnodes = linspace(t(1),t(end),N+1);
Bf = BernsteinPoly(f(tnodes),t);

tknots = [0 8 11 12];
N = 10;
M = length(tknots) - 1; 
tnodes = zeros(1,N*M+1);
for i = 1 : M
    tnodes(1,(i-1)*N+1:i*N+1) = linspace(tknots(i),tknots(i+1),N+1);
end
CoBf = CompositeBernsteinPoly(f(tnodes),tknots,t)


plot(t,f(t),'Color','b'); hold on
plot(t,Bf,'Color','r');
plot(t,CoBf,'Color','g');



% grid on

