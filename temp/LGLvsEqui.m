clear all
close all

N = 5; % order
tf = 3;



%% Compute Equidistant nodes from 0 to tf
t_equi = linspace(0,tf,N+1);

figure(1)
plot(t_equi,1,'o','Linewidth',3,'Color','r'); hold on
set(gca,'fontsize', 26);
grid on


%% Compute LGL nodes from 0 to tf
[t_LGL,w,D] = LGL_PS(N,tf);

figure(1)
plot(t_LGL,1,'o','Linewidth',3,'Color','g'); hold on
