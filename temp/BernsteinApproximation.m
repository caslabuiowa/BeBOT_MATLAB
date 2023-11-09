clear all 
close all

%% Define the following functions
T = 12;
u1 = @(t) 1./(sqrt(1+(T-t).^2));
u2 = @(t) (T-t)./(sqrt(1+(T-t).^2));
% Plot the functions in [0,T]
t = 0:0.01:T;
figure(1)
plot(t,u1(t),'LineWidth',3,'Color','k'); hold on
set(gca,'fontsize', 26);
grid on
figure(2)
plot(t,u2(t),'Linewidth',3,'Color','k'); hold on
grid on
set(gca,'fontsize', 26);




%% Approximation at equidistant nodes
N = 4;
[t_nodes, ~, ~] = BeBOT(N,T);
u1_nodes = u1(t_nodes);
u2_nodes = u2(t_nodes);

figure(1)
plot(t_nodes,u1_nodes,'o','Color','r','LineWidth',3);

figure(2)
plot(t_nodes,u2_nodes,'o','Color','r','LineWidth',3);

u1_N = BernsteinPoly(u1_nodes,t);
u2_N = BernsteinPoly(u2_nodes,t);

figure(1)
plot(t,u1_N,'Color','r','LineWidth',3);

figure(2)
plot(t,u2_N,'Color','r','LineWidth',3);



