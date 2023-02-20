clear all
close all

%% Load Par
CONSTANTS.N = 5; % Order of each poly
CONSTANTS.M = 2; % Number of polys
CONSTANTS.amax = 5;
CONSTANTS.amin = -5;
CONSTANTS.x10 = -3;
CONSTANTS.x1f = 0;
CONSTANTS.x20 = 0;
CONSTANTS.x2f = 0;




%% Init guess
N = CONSTANTS.N;
M = CONSTANTS.M; 
%x1 = linspace(CONSTANTS.x0,CONSTANTS.xf,CONSTANTS.N+1)';
x1 = ones((N+1)*M,1);
x2 = 2*ones((N+1)*M,1);
u = 3*ones((N+1)*M,1);
tf = 10;
x0 = [x1;x2;u;tf];


















%% Linear
A=[]; b=[]; Aeq=[]; beq=[]; lb=[]; ub=[]; 


















%% Optimize
options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',300000);
tic
[x,f] = fmincon(@(x)costfun(x,CONSTANTS),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,CONSTANTS),options);
toc






















%% Plot
N = CONSTANTS.N; 
M = CONSTANTS.M;
x1 = x(1:(N+1)*M);
x2 = x((N+1)*M+1:2*(N+1)*M);
u = x(2*(N+1)*M+1:3*(N+1)*M);
tf = x(end);
tknots = linspace(0,tf,M+1);
[tnodes,w,Diff] = PiecewiseBeBOT(N,tknots);
t = linspace(0,tf,1000);

figure
plot(tnodes,x1,'o'); hold on;
plot(tnodes,x2,'o');
plot(tknots(1:end-1),x1(1:N+1:end),'o','Linewidth',4,'Color','k');
plot(tknots(2:end),x1(N+1:N+1:end),'o','Linewidth',4,'Color','k');
plot(tknots(1:end-1),x2(1:N+1:end),'o','Linewidth',4,'Color','k');
plot(tknots(2:end),x2(N+1:N+1:end),'o','Linewidth',4,'Color','k');
plot(t,PiecewiseBernsteinPoly(x1',tknots,t));
plot(t,PiecewiseBernsteinPoly(x2',tknots,t));

figure
plot(tnodes,u,'o'); hold on;
plot(tknots(1:end-1),u(1:N+1:end),'o','Linewidth',4,'Color','k');
plot(tknots(2:end),u(N+1:N+1:end),'o','Linewidth',4,'Color','k');
plot(t,PiecewiseBernsteinPoly(u',tknots,t));



















%% Cost
function J = costfun(x,CONSTANTS)
%COSTFUN Summary of this function goes here
tf = x(end);

J = tf;
end














%% Constraints
function [c,ceq] = nonlcon(x,CONSTANTS)
%NONLCON Summary of this function goes here
%   Detailed explanation goes here
N = CONSTANTS.N; 
M = CONSTANTS.M;
x1 = x(1:(N+1)*M);
x2 = x((N+1)*M+1:2*(N+1)*M);
u = x(2*(N+1)*M+1:3*(N+1)*M);
tf = x(end);
tknots = linspace(0,tf,M+1);
[~,~,Diff] = PiecewiseBeBOT(CONSTANTS.N,tknots);

dyn1 = x1'*Diff-x2';
dyn2 = x2'*Diff-u';

nonlcon1 = u' - CONSTANTS.amax;
nonlcon2 = -u' + CONSTANTS.amin;

%% continuity constraints
cont_x1 = zeros(M-1,1);
cont_x2 = zeros(M-1,1);
if M > 1
    for i = 1 : M - 1
        cont_x1(i) = x1(i*(N+1)) - x1(i*(N+1)+1);  
        cont_x2(i) = x2(i*(N+1)) - x2(i*(N+1)+1);  
    end
end

c=[nonlcon1'; nonlcon2';];
ceq=[dyn1'; dyn2'; x1(1)-CONSTANTS.x10; x1(end)-CONSTANTS.x1f; ... 
        x2(1)-CONSTANTS.x20; x2(end)-CONSTANTS.x2f; ...
        cont_x1 ; cont_x2];
end