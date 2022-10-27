%% Bernstein Polynomials - Properties & Examples 
clear all
close all

tf = 10;
t = linspace(0,tf,1000);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples of 1,2,3D Bernstein polynomials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1D Bernstein Polynomial
x_N = [2 1 2 4 3];
xN = BernsteinPoly(x_N,t);


% 2D Bernstein Polynomials
x1_N = [6 0 5 9; 2 7 5 6];
x2_N = [0 -2 3 3; -4 1 -2 0];
x3_N = [6 4 10 9; 1 -2 -0 1];
x1N = BernsteinPoly(x1_N,t);
x2N = BernsteinPoly(x2_N,t);
x3N = BernsteinPoly(x3_N,t);


% 3D Bernstein Polynomials
p1_N = [1,2,4,5,7; 1,0,2,3,3; 1,1.2,2,3.5,3.5];
p2_N = [1,2,3,5,7; 7,7,6,4,4.5; 1,2,2.5,3,3];
p1N = BernsteinPoly(p1_N,t);
p2N = BernsteinPoly(p2_N,t);


% Plot the polynomials above
figure
subplot(2,2,1);
plot(t,xN,'LineWidth',3); hold on
plot(linspace(0,tf,length(x_N)),x_N,'o','LineWidth',3);
set(gca,'fontsize', 40);
axis([0 10 0 5])
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.52, 0.45, 0.45]);
legend('Bernstein polynomials','Bernstein coefficients')
subplot(2,2,3)
plot(x1N(1,:),x1N(2,:),'LineWidth',3); hold on
plot(x2N(1,:),x2N(2,:),'LineWidth',3);
plot(x3N(1,:),x3N(2,:),'LineWidth',3);
plot(x1_N(1,:),x1_N(2,:),'o','LineWidth',3);
plot(x2_N(1,:),x2_N(2,:),'o','LineWidth',3);
plot(x3_N(1,:),x3_N(2,:),'o','LineWidth',3);
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.02, 0.45, 0.45]);
subplot(2,2,[2,4])
plot3(p1N(1,:),p1N(2,:),p1N(3,:),'LineWidth',3); hold on
plot3(p2N(1,:),p2N(2,:),p2N(3,:),'LineWidth',3);
plot3(p1_N(1,:),p1_N(2,:),p1_N(3,:),'o','LineWidth',3);
plot3(p2_N(1,:),p2_N(2,:),p2_N(3,:),'o','LineWidth',3);
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.15, 0.8, 0.8]);






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depiction of End point value property
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot
figure
subplot(2,2,1);
plot(t,xN,'LineWidth',3); hold on
tnodes = linspace(0,tf,length(x_N));
plot(tnodes,x_N,'o','LineWidth',3);
plot(tnodes(1),x_N(1),'o','LineWidth',4,'Color','g');
plot(tnodes(end),x_N(end),'o','LineWidth',4,'Color','g');
plot(tnodes(1:2),x_N(1:2),'LineWidth',2.5,'Color','g');
plot(tnodes(end-1:end),x_N(end-1:end),'LineWidth',2.5,'Color','g');
set(gca,'fontsize', 40);
axis([0 10 0 5])
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.52, 0.45, 0.45]);
subplot(2,2,3)
plot(x1N(1,:),x1N(2,:),'LineWidth',3); hold on
plot(x2N(1,:),x2N(2,:),'LineWidth',3);
plot(x3N(1,:),x3N(2,:),'LineWidth',3);
plot(x1_N(1,:),x1_N(2,:),'o','LineWidth',3);
plot(x2_N(1,:),x2_N(2,:),'o','LineWidth',3);
plot(x3_N(1,:),x3_N(2,:),'o','LineWidth',3);
plot(x1_N(1,1),x1_N(2,1),'o','LineWidth',4,'Color','g');
plot(x1_N(1,end),x1_N(2,end),'o','LineWidth',4,'Color','g');
plot(x1_N(1,1:2),x1_N(2,1:2),'LineWidth',2.5,'Color','g');
plot(x1_N(1,end-1:end),x1_N(2,end-1:end),'LineWidth',2.5,'Color','g');
plot(x2_N(1,1),x2_N(2,1),'o','LineWidth',4,'Color','g');
plot(x2_N(1,end),x2_N(2,end),'o','LineWidth',4,'Color','g');
plot(x2_N(1,1:2),x2_N(2,1:2),'LineWidth',2.5,'Color','g');
plot(x2_N(1,end-1:end),x2_N(2,end-1:end),'LineWidth',2.5,'Color','g');
plot(x3_N(1,1),x3_N(2,1),'o','LineWidth',4,'Color','g');
plot(x3_N(1,end),x3_N(2,end),'o','LineWidth',4,'Color','g');
plot(x3_N(1,1:2),x3_N(2,1:2),'LineWidth',2.5,'Color','g');
plot(x3_N(1,end-1:end),x3_N(2,end-1:end),'LineWidth',2.5,'Color','g');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.02, 0.45, 0.45]);
subplot(2,2,[2,4])
plot3(p1N(1,:),p1N(2,:),p1N(3,:),'LineWidth',3); hold on
plot3(p2N(1,:),p2N(2,:),p2N(3,:),'LineWidth',3);
plot3(p1_N(1,:),p1_N(2,:),p1_N(3,:),'o','LineWidth',3);
plot3(p2_N(1,:),p2_N(2,:),p2_N(3,:),'o','LineWidth',3);
plot3(p1_N(1,1),p1_N(2,1),p1_N(3,1),'o','LineWidth',4,'Color','g');
plot3(p1_N(1,end),p1_N(2,end),p1_N(3,end),'o','LineWidth',4,'Color','g');
plot3(p1_N(1,1:2),p1_N(2,1:2),p1_N(3,1:2),'LineWidth',2.5,'Color','g');
plot3(p1_N(1,end-1:end),p1_N(2,end-1:end),p1_N(3,end-1:end),'LineWidth',2.5,'Color','g');
plot3(p2_N(1,1),p2_N(2,1),p2_N(3,1),'o','LineWidth',4,'Color','g');
plot3(p2_N(1,end),p2_N(2,end),p2_N(3,end),'o','LineWidth',4,'Color','g');
plot3(p2_N(1,1:2),p2_N(2,1:2),p2_N(3,1:2),'LineWidth',2.5,'Color','g');
plot3(p2_N(1,end-1:end),p2_N(2,end-1:end),p2_N(3,end-1:end),'LineWidth',2.5,'Color','g');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.15, 0.8, 0.8]);





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Depiction of convex hull property
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot
figure
subplot(2,2,1);
plot(t,xN,'LineWidth',3); hold on
plot(linspace(0,tf,length(x_N)),x_N,'o','LineWidth',3);
plot(linspace(0,tf,length(x_N)),ones(length(x_N))*max(x_N),'--','LineWidth',0.5,'Color','r');
plot(linspace(0,tf,length(x_N)),ones(length(x_N))*min(x_N),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis([0 10 0 5])
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.52, 0.45, 0.45]);
subplot(2,2,3)
plot(x1N(1,:),x1N(2,:),'LineWidth',3); hold on
plot(x2N(1,:),x2N(2,:),'LineWidth',3);
plot(x3N(1,:),x3N(2,:),'LineWidth',3);
plot(x1_N(1,:),x1_N(2,:),'o','LineWidth',3);
plot(x2_N(1,:),x2_N(2,:),'o','LineWidth',3);
plot(x3_N(1,:),x3_N(2,:),'o','LineWidth',3);
[k,av] = convhull(x1_N');
plot(x1_N(1,k),x1_N(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x2_N');
plot(x2_N(1,k),x2_N(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x3_N');
plot(x3_N(1,k),x3_N(2,k),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.02, 0.45, 0.45]);
subplot(2,2,[2,4])
plot3(p1N(1,:),p1N(2,:),p1N(3,:),'LineWidth',3); hold on
plot3(p2N(1,:),p2N(2,:),p2N(3,:),'LineWidth',3);
plot3(p1_N(1,:),p1_N(2,:),p1_N(3,:),'o','LineWidth',3);
plot3(p2_N(1,:),p2_N(2,:),p2_N(3,:),'o','LineWidth',3);
[k,av] = convhull(p1_N');
trisurf(k,p1_N(1,:),p1_N(2,:),p1_N(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
[k,av] = convhull(p2_N');
trisurf(k,p2_N(1,:),p2_N(2,:),p2_N(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.15, 0.8, 0.8]);








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Degree evelation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = 10; % Elevate to Mth order

% Degree elevate the coefficients of the 1D Bernstein polynomial
x_N_elev = x_N*DegElevMatrix(length(x_N)-1,M);

% Degree elevate the coefficients of the 2D Bernstein polynomials
x1_N_elev = x1_N*DegElevMatrix(length(x1_N)-1,M);
x2_N_elev = x2_N*DegElevMatrix(length(x2_N)-1,M);
x3_N_elev = x3_N*DegElevMatrix(length(x3_N)-1,M);

% Degree elevate the coefficients of the 3D Bernstein polynomials
p1_N_elev = p1_N*DegElevMatrix(length(p1_N)-1,M);
p2_N_elev = p2_N*DegElevMatrix(length(p2_N)-1,M);


% Plot
figure
subplot(2,2,1);
plot(t,xN,'LineWidth',3); hold on
plot(linspace(0,tf,length(x_N_elev)),x_N_elev,'o','LineWidth',3);
plot(linspace(0,tf,length(x_N_elev)),ones(length(x_N_elev))*max(x_N_elev),'--','LineWidth',0.5,'Color','r');
plot(linspace(0,tf,length(x_N_elev)),ones(length(x_N_elev))*min(x_N_elev),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis([0 10 0 5])
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.52, 0.45, 0.45]);
subplot(2,2,3)
plot(x1N(1,:),x1N(2,:),'LineWidth',3); hold on
plot(x2N(1,:),x2N(2,:),'LineWidth',3);
plot(x3N(1,:),x3N(2,:),'LineWidth',3);
plot(x1_N_elev(1,:),x1_N_elev(2,:),'o','LineWidth',3);
plot(x2_N_elev(1,:),x2_N_elev(2,:),'o','LineWidth',3);
plot(x3_N_elev(1,:),x3_N_elev(2,:),'o','LineWidth',3);
[k,av] = convhull(x1_N_elev');
plot(x1_N_elev(1,k),x1_N_elev(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x2_N_elev');
plot(x2_N_elev(1,k),x2_N_elev(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x3_N_elev');
plot(x3_N_elev(1,k),x3_N_elev(2,k),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.02, 0.45, 0.45]);
subplot(2,2,[2,4])
plot3(p1N(1,:),p1N(2,:),p1N(3,:),'LineWidth',3); hold on
plot3(p2N(1,:),p2N(2,:),p2N(3,:),'LineWidth',3);
plot3(p1_N_elev(1,:),p1_N_elev(2,:),p1_N_elev(3,:),'o','LineWidth',3);
plot3(p2_N_elev(1,:),p2_N_elev(2,:),p2_N_elev(3,:),'o','LineWidth',3);
[k,av] = convhull(p1_N_elev');
trisurf(k,p1_N_elev(1,:),p1_N_elev(2,:),p1_N_elev(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
[k,av] = convhull(p2_N_elev');
trisurf(k,p2_N_elev(1,:),p2_N_elev(2,:),p2_N_elev(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.15, 0.8, 0.8]);









%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% de Casteljau algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 0.5; % Value at which we subdivide (lambda \in [0,1])

% Subdivide the 1D Bernstein polynomials at lambda 
[Cpout Pos] = deCasteljau(x_N,lambda); 
% Pos is the value of xN at lambda*tf (or xN(lambda*tf))
x_NA = Cpout(:, 1:length(x_N(1,:)));
x_NB = Cpout(:, length(x_N(1,:)):end);
% x_NA and x_NB are the coefficients of the two Bernstein polynomials
xNA = BernsteinPoly(x_NA,t);
xNB = BernsteinPoly(x_NB,t);

% Plot
figure
subplot(2,2,1);
tA = linspace(0,lambda,length(xNA))*tf;
tB = linspace(lambda,1,length(xNA))*tf;
plot(tA,xNA,'LineWidth',3); hold on
plot(tB,xNB,'LineWidth',3);
plot(linspace(0,lambda*tf,length(x_NA)),x_NA,'o','LineWidth',3);
plot(linspace(lambda*tf,tf,length(x_NB)),x_NB,'o','LineWidth',3);
plot(lambda*tf,Pos,'o','LineWidth',5,'Color','r');
plot(linspace(0,tf,length(x_N)),ones(length(x_N))*max(max(x_NA),max(x_NB)),'--','LineWidth',0.5,'Color','r');
plot(linspace(0,tf,length(x_N)),ones(length(x_N))*min(min(x_NA),min(x_NB)),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis([0 10 0 5])
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.52, 0.45, 0.45]);



% Subdivide the 2D Bernstein polynomials at lambda 
[Cpout Pos1] = deCasteljau(x1_N,lambda);
x1_NA = Cpout(:, 1:length(x1_N(1,:)));
x1_NB = Cpout(:, length(x1_N(1,:)):end);
x1NA = BernsteinPoly(x1_NA,t);
x1NB = BernsteinPoly(x1_NB,t);
[Cpout Pos2] = deCasteljau(x2_N,lambda);
x2_NA = Cpout(:, 1:length(x2_N(1,:)));
x2_NB = Cpout(:, length(x2_N(1,:)):end);
x2NA = BernsteinPoly(x2_NA,t);
x2NB = BernsteinPoly(x2_NB,t);
[Cpout Pos3] = deCasteljau(x3_N,lambda);
x3_NA = Cpout(:, 1:length(x3_N(1,:)));
x3_NB = Cpout(:, length(x3_N(1,:)):end);
x3NA = BernsteinPoly(x3_NA,t);
x3NB = BernsteinPoly(x3_NB,t);

% Plot
subplot(2,2,3)
plot(x1NA(1,:),x1NA(2,:),'LineWidth',3); hold on
plot(x1NB(1,:),x1NB(2,:),'LineWidth',3); 
plot(x2NA(1,:),x2NA(2,:),'LineWidth',3);
plot(x2NB(1,:),x2NB(2,:),'LineWidth',3);
plot(x3NA(1,:),x3NA(2,:),'LineWidth',3);
plot(x3NB(1,:),x3NB(2,:),'LineWidth',3);
plot(x1_NA(1,:),x1_NA(2,:),'o','LineWidth',3);
plot(x2_NA(1,:),x2_NA(2,:),'o','LineWidth',3);
plot(x3_NA(1,:),x3_NA(2,:),'o','LineWidth',3);
plot(x1_NB(1,:),x1_NB(2,:),'o','LineWidth',3);
plot(x2_NB(1,:),x2_NB(2,:),'o','LineWidth',3);
plot(x3_NB(1,:),x3_NB(2,:),'o','LineWidth',3);
plot(Pos1(1),Pos1(2),'o','LineWidth',5,'Color','r');
plot(Pos2(1),Pos2(2),'o','LineWidth',5,'Color','r');
plot(Pos3(1),Pos3(2),'o','LineWidth',5,'Color','r');
[k,av] = convhull(x1_NA');
plot(x1_NA(1,k),x1_NA(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x1_NB');
plot(x1_NB(1,k),x1_NB(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x2_NA');
plot(x2_NA(1,k),x2_NA(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x2_NB');
plot(x2_NB(1,k),x2_NB(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x3_NA');
plot(x3_NA(1,k),x3_NA(2,k),'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(x3_NB');
plot(x3_NB(1,k),x3_NB(2,k),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.02, 0.45, 0.45]);



% Subdivide the 3D Bernstein polynomials at lambda
[Cpout Pos1] = deCasteljau(p1_N,lambda);
p1_NA = Cpout(:, 1:length(p1_N(1,:)));
p1_NB = Cpout(:, length(p1_N(1,:)):end);
p1NA = BernsteinPoly(p1_NA,t);
p1NB = BernsteinPoly(p1_NB,t);
[Cpout Pos2] = deCasteljau(p2_N,lambda);
p2_NA = Cpout(:, 1:length(p2_N(1,:)));
p2_NB = Cpout(:, length(p2_N(1,:)):end);
p2NA = BernsteinPoly(p2_NA,t);
p2NB = BernsteinPoly(p2_NB,t);
subplot(2,2,[2,4])
plot3(p1NA(1,:),p1NA(2,:),p1NA(3,:),'LineWidth',3); hold on
plot3(p1NB(1,:),p1NB(2,:),p1NB(3,:),'LineWidth',3);
plot3(p2NA(1,:),p2NA(2,:),p2NA(3,:),'LineWidth',3);
plot3(p2NB(1,:),p2NB(2,:),p2NB(3,:),'LineWidth',3);
plot3(p1_NA(1,:),p1_NA(2,:),p1_NA(3,:),'o','LineWidth',3);
plot3(p1_NB(1,:),p1_NB(2,:),p1_NB(3,:),'o','LineWidth',3);
plot3(p2_NA(1,:),p2_NA(2,:),p2_NA(3,:),'o','LineWidth',3);
plot3(p2_NB(1,:),p2_NB(2,:),p2_NB(3,:),'o','LineWidth',3);
plot3(Pos1(1),Pos1(2),Pos1(3),'o','LineWidth',5,'Color','r');
plot3(Pos2(1),Pos2(2),Pos2(3),'o','LineWidth',5,'Color','r');
[k,av] = convhull(p1_NA');
trisurf(k,p1_NA(1,:),p1_NA(2,:),p1_NA(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
[k,av] = convhull(p1_NB');
trisurf(k,p1_NB(1,:),p1_NB(2,:),p1_NB(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
[k,av] = convhull(p2_NA');
trisurf(k,p2_NA(1,:),p2_NA(2,:),p2_NA(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
[k,av] = convhull(p2_NB');
trisurf(k,p2_NB(1,:),p2_NB(2,:),p2_NB(3,:),'FaceColor','r','FaceAlpha',.1,'EdgeAlpha',.1)
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.15, 0.8, 0.8]);







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Minimum distance algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute max and min of a 1D Bernstein polynomial
[max, tmax] = MaximumBernstein(x_N);
[~, pmax] = deCasteljau(x_N,tmax);
[min, tmin] = MinimumBernstein(x_N);
[~, pmin] = deCasteljau(x_N,tmin);

% Plot
figure
subplot(2,2,1);
plot(t,xN,'LineWidth',3); hold on
tnodes = linspace(0,tf,length(x_N));
plot(tmax*tf,max,'o','LineWidth',4,'Color','g');
plot(linspace(0,10,tf),max*ones(10,1),'--','LineWidth',0.5,'Color','r');
plot(tmin*tf,min,'o','LineWidth',4,'Color','g');
plot(linspace(0,10,tf),min*ones(10,1),'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis([0 10 0 5])
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.52, 0.45, 0.45]);

% Compute distance between 2D Bernstein polynomials
% and between 2D Bernstein polynomials and convex shapes
[dist, t1, t2] = MinDistBernstein2Bernstein(x1_N, x2_N);
[~,x12col1] = deCasteljau(x1_N,t1);
[~,x12col2] = deCasteljau(x2_N,t2);
[dist, t1, t3] = MinDistBernstein2Bernstein(x1_N, x3_N);
[~,x13col1] = deCasteljau(x1_N,t1);
[~,x13col3] = deCasteljau(x3_N,t3);
[dist, t2, t3] = MinDistBernstein2Bernstein(x2_N, x3_N);
[~,x23col2] = deCasteljau(x2_N,t2);
[~,x23col3] = deCasteljau(x3_N,t3);
shape = [1 -1  -0.2; 4 6 2];
[dist, t, pt] = MinDistBernstein2Polygon(x1_N, shape);
[~,x1shapecol1] = deCasteljau(x1_N,t);

% Plot
subplot(2,2,3)
plot(x1N(1,:),x1N(2,:),'LineWidth',3); hold on
plot(x2N(1,:),x2N(2,:),'LineWidth',3);
plot(x3N(1,:),x3N(2,:),'LineWidth',3);
plot(x12col1(1),x12col1(2),'o','LineWidth',4,'Color','r');
plot(x12col2(1),x12col2(2),'o','LineWidth',4,'Color','r');
plot([x12col1(1) x12col2(1)],[x12col1(2) x12col2(2)],'--','LineWidth',0.5,'Color','r');
plot(x13col1(1),x13col1(2),'o','LineWidth',4,'Color','r');
plot(x13col3(1),x13col3(2),'o','LineWidth',4,'Color','r');
plot([x13col1(1) x13col3(1)],[x13col1(2) x13col3(2)],'--','LineWidth',0.5,'Color','r');
plot(x23col2(1),x23col2(2),'o','LineWidth',4,'Color','r');
plot(x23col3(1),x23col3(2),'o','LineWidth',4,'Color','r');
plot([x23col2(1) x23col3(1)],[x23col2(2) x23col3(2)],'--','LineWidth',0.5,'Color','r');
plot(x1shapecol1(1),x1shapecol1(2),'o','LineWidth',4,'Color','r');
plot(pt(1),pt(2),'o','LineWidth',4,'Color','r');
plot([x1shapecol1(1) pt(1)],[x1shapecol1(2) pt(2)],'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(shape');
plot(shape(1,k),shape(2,k),'-','LineWidth',5,'Color','k');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.02, 0.02, 0.45, 0.45]);

% Compute distance between 3D Bernstein polynomials
[dist, t1, t2] = MinDistBernstein2Bernstein(p1_N, p2_N);
[~,pcol1] = deCasteljau(p1_N,t1);
[~,pcol2] = deCasteljau(p2_N,t2);
shape = [3 3 3 5; 4 4 5 6; 0 2 2 1.5];
[dist, t, pt1] = MinDistBernstein2Polygon(p1_N, shape);
[~,pcol1shape] = deCasteljau(p1_N,t);
[dist, t, pt2] = MinDistBernstein2Polygon(p2_N, shape);
[~,pcol2shape] = deCasteljau(p2_N,t);

% Plot
subplot(2,2,[2,4])
plot3(p1N(1,:),p1N(2,:),p1N(3,:),'LineWidth',3); hold on
plot3(p2N(1,:),p2N(2,:),p2N(3,:),'LineWidth',3);
plot3(pcol1(1),pcol1(2),pcol1(3),'o','LineWidth',4,'Color','r');
plot3(pcol2(1),pcol2(2),pcol2(3),'o','LineWidth',4,'Color','r');
plot3([pcol1(1) pcol2(1)],[pcol1(2) pcol2(2)],[pcol1(3) pcol2(3)],'--','LineWidth',0.5,'Color','r');
set(gca,'fontsize', 40);
axis equal
grid on
set(gca, 'Units', 'Normalized', 'OuterPosition', [0.35, 0.15, 0.8, 0.8]);
plot3(pcol1shape(1),pcol1shape(2),pcol1shape(3),'o','LineWidth',4,'Color','r');
plot3(pcol2shape(1),pcol2shape(2),pcol2shape(3),'o','LineWidth',4,'Color','r');
plot3(pt1(1),pt1(2),pt1(3),'o','LineWidth',4,'Color','r');
plot3(pt2(1),pt2(2),pt2(3),'o','LineWidth',4,'Color','r');
plot3([pcol1shape(1) pt1(1)],[pcol1shape(2) pt1(2)],[pcol1shape(3) pt1(3)],'--','LineWidth',0.5,'Color','r');
plot3([pcol2shape(1) pt2(1)],[pcol2shape(2) pt2(2)],[pcol2shape(3) pt2(3)],'--','LineWidth',0.5,'Color','r');
[k,av] = convhull(shape');
trisurf(k,shape(1,:),shape(2,:),shape(3,:),'FaceColor','k','FaceAlpha',.5,'EdgeAlpha',.5)
set(gca,'fontsize', 40);


%% CheckColDist: check if two Bernstein polynomials are further than mindist
% apart
collision = CollCheckBernstein2Bernstein(p1_N, p2_N, 1.6)
collision = CollCheckBernstein2Bernstein(p1_N, p2_N, 1.5)
collision = CollCheckBernstein2Polygon(p1_N, shape, 2.4)
collision = CollCheckBernstein2Polygon(p1_N, shape, 2.3)
%[dist, t1, t2] = MinDistBernstein2Bernstein(p1_N, p2_N)






























