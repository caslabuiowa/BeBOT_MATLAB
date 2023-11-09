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


% Their integrals are
intu1 = asinh(T);
intu2 = sqrt(T^2+1)-1;

%% Numerical Integration
% Compute weights
N = 70;
[tnodes, w, D] = BeBOT(N,T);
% Evaluate function at the nodes
u1_N = u1(tnodes);
figure(1) 
plot(tnodes,u1_N,'o','LineWidth',3,'Color','k'); hold on
u2_N = u2(tnodes);
figure(2) 
plot(tnodes,u2_N,'o','LineWidth',3,'Color','k'); hold on

int_u1 = u1_N*w;
int_u2 = u2_N*w;

err1 = int_u1 - intu1
err2 = int_u2 - intu2



