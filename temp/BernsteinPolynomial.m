clear all
close all

time = 0:0.01:7;

xbar = [0 -2  -1  3];

t_nodes = linspace(0,7,length(xbar));

xN = BernsteinPoly(xbar,time);

figure(1)
plot(t_nodes,xbar,'o','Linewidth',3,'Color','b'); hold on
plot(time,xN,'Linewidth',3,'Color','r');



%%
for i = 0 : 500

    xbar = [0 -2+0.01*i  -1  3];

    t_nodes = linspace(0,7,length(xbar));

    xN = BernsteinPoly(xbar,time);

    figure(2)
    hold off
    plot(t_nodes,xbar,'o','Linewidth',3,'Color','b'); hold on
    plot(time,xN,'Linewidth',3,'Color','r');
    axis([0 7 -2 3])

    pause(0.01)

end

for i = 0 : 400

    xbar = [0 -2+5  -1+0.01*i  3];

    t_nodes = linspace(0,7,length(xbar));

    xN = BernsteinPoly(xbar,time);

    figure(2)
    hold off
    plot(t_nodes,xbar,'o','Linewidth',3,'Color','b'); hold on
    plot(time,xN,'Linewidth',3,'Color','r');
    axis([0 7 -2 3])

    pause(0.01)

end

%% Higher order 

xbar = [0 -0.4 -0.8 -1.2 -1.4 -1.7 -2 -1.8 -1.6 -1.4 -1.2 -1 -0.5 0 0.5 1 1.5 2 2.5 3];

t_nodes = linspace(0,7,length(xbar));

xN = BernsteinPoly(xbar,time);

figure(3)
hold off
plot(t_nodes,xbar,'o','Linewidth',3,'Color','b'); hold on
plot(time,xN,'Linewidth',3,'Color','r');
axis([0 7 -2 3])

%% Bernstein and Lagrange
xbar = [0 -2  -1  3];
t_nodes = linspace(0,7,length(xbar));
xN = BernsteinPoly(xbar,time);
figure(4)
plot(t_nodes,xbar,'o','Linewidth',3,'Color','b'); hold on
plot(time,xN,'Linewidth',3,'Color','r');
xN = LagrangePoly(xbar,t_nodes,time);
plot(time,xN,'Linewidth',3,'Color','g');
