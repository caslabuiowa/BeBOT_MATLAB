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

% Their derivatives are
u1dot = @(t) (T-t)./((T-t).^2+1).^(3/2);
u2dot = @(t) -1./(T^2-2*T*t+t.^2+1).^(3/2);
% Plot the derivatives in [0,T]
figure(3)
plot(t,u1dot(t),'Linewidth',3,'Color','k'); hold on
set(gca,'fontsize', 26);
grid on
figure(4)
plot(t,u2dot(t),'Linewidth',3,'Color','k'); hold on
set(gca,'fontsize', 26);
grid on





%% Approximation at LGL nodes
N = 10;
[t_nodes,w,Dm] = LGL_PS(N,T);
u1_nodes = u1(t_nodes);
u2_nodes = u2(t_nodes);

figure(1)
plot(t_nodes,u1_nodes,'o','Color','g','LineWidth',3);

figure(2)
plot(t_nodes,u2_nodes,'o','Color','g','LineWidth',3);

u1_N = LagrangePoly(u1_nodes,t_nodes,t);
u2_N = LagrangePoly(u2_nodes,t_nodes,t);

figure(1)
plot(t,u1_N,'Color','g','LineWidth',3);

figure(2)
plot(t,u2_N,'Color','g','LineWidth',3);



%% Approximate the derivative

u1dot_nodes = u1_nodes*Dm;
u2dot_nodes = u2_nodes*Dm;

figure(3)
plot(t_nodes,u1dot_nodes,'o','Color','g','LineWidth',3);

figure(4)
plot(t_nodes,u2dot_nodes,'o','Color','g','LineWidth',3);

u1dot_N = LagrangePoly(u1dot_nodes,t_nodes,t);
u2dot_N = LagrangePoly(u2dot_nodes,t_nodes,t);

figure(3)
plot(t,u1dot_N,'Color','g','LineWidth',3);

figure(4)
plot(t,u2dot_N,'Color','g','LineWidth',3);

