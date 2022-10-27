clear all
close all

time = 0:0.01:7;

t_nodes = [0 1 3 5 7];
xbar = [2 1 -2 -1 3];

xN = LagrangePoly(xbar,t_nodes,time);

figure(1)
plot(t_nodes,xbar,'o','Linewidth',3,'Color','b'); hold on
plot(time,xN,'Linewidth',3,'Color','r');


pbar = [2 1 -2 -1 3; ...
    2 5 9 20 25];

pNx = LagrangePoly(pbar(1,:),t_nodes,time);
pNy = LagrangePoly(pbar(2,:),t_nodes,time);
figure(2)
plot(pbar(1,:),pbar(2,:),'o','Linewidth',3,'Color','b'); hold on
plot(pNx,pNy,'Linewidth',3,'Color','r');