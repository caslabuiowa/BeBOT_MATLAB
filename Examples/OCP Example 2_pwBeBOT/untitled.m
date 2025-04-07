clear all
close all


CONSTANTS.N = 3; % Order of approximation
CONSTANTS.M = 2; % Number of knots
M = CONSTANTS.M;
CONSTANTS.tknots = linspace(0,1,M);
CONSTANTS.amax = 5;
CONSTANTS.amin = -5;
CONSTANTS.x10 = -3;
CONSTANTS.x1f = 0;
CONSTANTS.x20 = 0;
CONSTANTS.x2f = 0;



N = CONSTANTS.N; 
M = CONSTANTS.M; 
x1 = ones((N+1)*(M-1),1);
x2 = ones((N+1)*(M-1),1);
u = ones((N+1)*(M-1),1);
tf = 10;
x0 = [x1;x2;u;tf];



A=[]; b=[]; lb=[]; ub=[];

% Initialize Aeq and beq
Aeq = zeros(numel(x0), numel(x0));  % Create a zero matrix of appropriate size
beq = zeros(numel(x0), 1);  % Zero vector for equality constraints

if CONSTANTS.M > 2
    % Loop over each internal knot to enforce continuity
    for i = 1:(CONSTANTS.M - 2)  % M-1 internal knots
        % Indices for the last coefficient of x1 at the i-th segment
        idx_last_x1 = (N+1) * i;
        % Indices for the first coefficient of x1 at the (i+1)-th segment
        idx_first_x1_next = (N+1) * i + 1;

        % Indices for the last coefficient of x2 at the i-th segment
        idx_last_x2 = (N+1) * i + (N+1) * (M-1);
        % Indices for the first coefficient of x2 at the (i+1)-th segment
        idx_first_x2_next = (N+1) * i + (N+1) * (M-1) + 1;

        % Enforcing continuity of x1
        Aeq(idx_last_x1, idx_last_x1) = 1;
        Aeq(idx_last_x1, idx_first_x1_next) = -1;
        
        % Enforcing continuity of x2
        Aeq(idx_last_x2, idx_last_x2) = 1;
        Aeq(idx_last_x2, idx_first_x2_next) = -1;
    end
end



options = optimoptions(@fmincon,'Algorithm','sqp','MaxFunctionEvaluations',300000);
tic
[x,f] = fmincon(@(x)costfun(x,CONSTANTS),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x,CONSTANTS),options);
toc


N = CONSTANTS.N;
M = CONSTANTS.M;

x1 = x(1:(N+1)*(M-1));
x2 = x((N+1)*(M-1) + 1:2*(N+1)*(M-1));
u = x(2*(N+1)*(M-1) + 1 :3*(N+1)*(M-1));
tf = x(end);

tknots = CONSTANTS.tknots*tf;
[tnodes,w,Diff] = PiecewiseBeBOT(N,tknots);
t = linspace(0,tf,1000);

figure
plot(tnodes,x1,'o'); hold on
plot(tnodes,x2,'o');
plot(t,PiecewiseBernsteinPoly(x1',tknots,t));
plot(t,PiecewiseBernsteinPoly(x2',tknots,t));

figure
plot(tnodes,u,'o'); hold on
plot(t,PiecewiseBernsteinPoly(u',tknots,t));



function J = costfun(x,CONSTANTS)
%COSTFUN Summary of this function goes here
N = CONSTANTS.N;
M = CONSTANTS.M;

x1 = x(1:(N+1)*(M-1));
x2 = x((N+1)*(M-1) + 1:2*(N+1)*(M-1));
u = x(2*(N+1)*(M-1) + 1 :3*(N+1)*(M-1));
tf = x(end);

J = tf;
end


function [c,ceq] = nonlcon(x,CONSTANTS)
%NONLCON Summary of this function goes here
%   Detailed explanation goes here
N = CONSTANTS.N;
M = CONSTANTS.M;

x1 = x(1:(N+1)*(M-1));
x2 = x((N+1)*(M-1) + 1:2*(N+1)*(M-1));
u = x(2*(N+1)*(M-1) + 1 :3*(N+1)*(M-1));
tf = x(end);

tknots = CONSTANTS.tknots*tf;
[~,~,Diff] = PiecewiseBeBOT(CONSTANTS.N,tknots);

dyn1 = x1'*Diff-x2';
dyn2 = x2'*Diff-u';

nonlcon1 = u' - CONSTANTS.amax;
nonlcon2 = -u' + CONSTANTS.amin;

c=[nonlcon1'; nonlcon2';];
ceq=[dyn1'; dyn2'; x1(1)-CONSTANTS.x10; x1(end)-CONSTANTS.x1f; x2(1)-CONSTANTS.x20; x2(end)-CONSTANTS.x2f];
end