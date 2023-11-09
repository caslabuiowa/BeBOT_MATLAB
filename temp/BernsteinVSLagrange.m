clear all 
close all

N = 50;

%% Define the following functions
T = 2*pi;
f = @(t) cos(t);
t = 0:0.01:T;
figure(1)
plot(t,f(t),'LineWidth',3,'Color','k'); hold on
figure(2)
plot(t,f(t),'LineWidth',3,'Color','k'); hold on

%% Bernstein approximation at equidistant nodes
[t_nodes, ~, ~] = BeBOT(N,T);
f_nodes = f(t_nodes);

figure(1)
plot(t_nodes,f_nodes,'o','Color','r','LineWidth',3);

f_N = BernsteinPoly(f_nodes,t);

figure(1)
plot(t,f_N,'Color','r','LineWidth',3);


%% Lagrange approximation at LGL nodes
[t_nodes, ~, ~] = LGL_PS(N,T);
f_nodes = f(t_nodes);

figure(2)
plot(t_nodes,f_nodes,'o','Color','r','LineWidth',3);

f_N = LagrangePoly(f_nodes,t_nodes,t);

figure(2)
plot(t,f_N,'Color','r','LineWidth',3);
