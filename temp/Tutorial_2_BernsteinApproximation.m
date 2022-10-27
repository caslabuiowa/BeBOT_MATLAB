clear all
close all






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernstein Approximation of smooth functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 12; % Order of approximation (N+1 nodes)

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

% Their integrals are
intu1 = asinh(T);
intu2 = sqrt(T^2+1)-1;




%% Compute Bern approx nodes/weights/Diff
[tnodes,w,D] = BeBOT(N,T)


%% Plot the function at the equidistant nodes
u1nodes = u1(tnodes);
u2nodes = u2(tnodes);
figure(1)
plot(tnodes,u1nodes,'o','Linewidth',4); 
figure(2)
plot(tnodes,u2nodes,'o','Linewidth',4); 


%% Compute the Bernstein Polynomial at these nodes
u1N = BernsteinPoly(u1nodes,t);
u2N = BernsteinPoly(u2nodes,t);
figure(1)
plot(t,u1N,'Linewidth',2,'Color','b'); 
figure(2)
plot(t,u2N,'Linewidth',2,'Color','b'); 


%% Compute the derivative of Bernstein poly
u1nodes_dot = u1nodes*D;
u2nodes_dot = u2nodes*D;
u1Ndot = BernsteinPoly(u1nodes_dot,t);
u2Ndot = BernsteinPoly(u2nodes_dot,t);
figure(3)
plot(tnodes,u1nodes_dot,'o','Linewidth',4);
plot(t,u1Ndot,'Linewidth',2,'Color','b'); 
figure(4)
plot(tnodes,u2nodes_dot,'o','Linewidth',4); 
plot(t,u2Ndot,'Linewidth',2,'Color','b'); 


%% Compute the integral of Bernstein poly
intu1N = u1nodes*w;
intu2N = u2nodes*w;
% Compare to analytical answer
intu1;
intu2;


% Compare the integral of Bernstein approximation for different N
i = 0;
for N = 5:200
    i = i + 1;
    [tnodes,w,~] = BeBOT(N,T)
    u1nodes = u1(tnodes);
    u2nodes = u2(tnodes);
    intu1N(i) = u1nodes*w;
    intu2N(i) = u2nodes*w;
    order(i) = N;
end
figure(5)
plot(order,intu1N,'Linewidth',2,'Color','b'); hold on
plot(order,intu1*ones(length(order),1),'Linewidth',2,'Color','k');
legend('numerical','analytical')
set(gca,'fontsize', 26);
grid on
figure(6)
plot(order,intu2N,'Linewidth',2,'Color','b'); hold on
plot(order,intu2*ones(length(order),1),'Linewidth',2,'Color','k');
legend('numerical','analytical')
set(gca,'fontsize', 26);
grid on









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bernstein Approximation of non-smooth functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 15; % Order of approximation (N+1 nodes)


%% Define the following function
syms u(t)
T = 2;
u(t) =@(t) piecewise(t<=T/2, 1, T/2<t, -1);
t = 0:0.01:T;
uval = u(t);
figure(7)
plot(t,uval,'Linewidth',2,'Color','k'); hold on
set(gca,'fontsize', 26);
grid on

%% Compute Bern approx nodes/weights/Diff
[tnodes,~,~] = BeBOT(N,T);

%% Plot the function at the equidistant nodes
figure(7)
unodes = u(tnodes);
plot(tnodes,unodes,'o','Linewidth',4);


%% Compute the Bernstein Polynomial at these nodes
uN = BernsteinPoly(unodes,t);
figure(7)
plot(t,uN,'Linewidth',2,'Color','b'); 









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lagrange approximation of non-smooth functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the function on different figures
figure(8)
plot(t,uval,'Linewidth',2,'Color','k'); hold on
set(gca,'fontsize', 26);
grid on

%% Compute LGL nodes
[tnodes,~,~] = LGL_PS(N,T); % Compute tnodes,w,Diff in [0,T]

%% Plot the function at LGL nodes
figure(8)
unodes = u(tnodes);
plot(tnodes,unodes,'o','Linewidth',4);

%% Compute the Lagrange polynomial at these nodes
uN = LagrangePoly(unodes,tnodes,t);
figure(8)
plot(t,uN,'Linewidth',2,'Color','b'); 

