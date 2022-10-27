clear all
close all


%% Define the following function
syms u(t)
T = 2;
u(t) =@(t) piecewise(t<=T/2, 1, T/2<t, -1);
t = 0:0.01:T;
uval = u(t);
figure(1)
plot(t,uval,'Linewidth',2,'Color','k'); hold on
set(gca,'fontsize', 26);
grid on


%% Approximation at LGL nodes
N = 100;
[t_nodes,w,Dm] = LGL_PS(N,T)
u_nodes = u(t_nodes);


figure(1)
plot(t_nodes,u_nodes,'o','Color','g','LineWidth',3);


u_N = LagrangePoly(u_nodes,t_nodes,t);

figure(1)
plot(t,u_N,'Color','g','LineWidth',3);
