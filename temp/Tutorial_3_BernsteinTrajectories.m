clear all
close all

tf = 10;
t = linspace(0,tf,1000);



%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Motion Planning tools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2D Bernstein Polynomials
x1_N = [6 0 5 9; 2 7 5 6];
x2_N = [0 -2 3 3; -4 1 -2 0];
x1N = BernsteinPoly(x1_N,t);
x2N = BernsteinPoly(x2_N,t);


% Plot
figure(1)
plot(x1N(1,:),x1N(2,:),'LineWidth',3); hold on
plot(x2N(1,:),x2N(2,:),'LineWidth',3);
plot(x1_N(1,:),x1_N(2,:),'o','LineWidth',3);
plot(x2_N(1,:),x2_N(2,:),'o','LineWidth',3);
set(gca,'fontsize', 40);
axis equal
grid on


%% Speed square
[~,N] = size(x2_N);
N = N - 1;
[tnodes,w,Dm] = BeBOT(N,tf);
x2_Ndot = x2_N*Dm;
speedsqr_2N = BernsteinProduct(x2_Ndot(1,:),x2_Ndot(1,:))+BernsteinProduct(x2_Ndot(2,:),x2_Ndot(2,:));
speedsqr2N = BernsteinPoly(speedsqr_2N,t);
speed2N = sqrt(speedsqr2N);


% Plot
figure(2)
plot(t,speed2N,'LineWidth',3); hold on
plot(t,speedsqr2N,'LineWidth',3);
plot(linspace(0,tf,2*N+1),speedsqr_2N,'o','LineWidth',3);
set(gca,'fontsize', 40);
legend('speed','speedsqr','speedsqr coef')
grid on



%% Psidot
x2_Nddot = x2_Ndot*Dm;
psidot_num = (BernsteinProduct(x2_Nddot(1,:),x2_Ndot(2,:))-BernsteinProduct(x2_Nddot(2,:),x2_Ndot(1,:)));
psidot_num = psidot_num*DegElevMatrix(length(psidot_num)-1,length(psidot_num));
psidot_den = speedsqr_2N*DegElevMatrix(length(speedsqr_2N)-1,length(speedsqr_2N));
psidot_coeff = psidot_num./(psidot_den);
psidot_weights = (psidot_den); 
psidot = BernsteinPoly(psidot_num,t)./BernsteinPoly(psidot_den,t);



% Plot
figure(3)
plot(t,psidot,'LineWidth',3); hold on
plot(linspace(0,tf,length(psidot_coeff)),psidot_coeff,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on




% Degree elevate
[~,N] = size(psidot_num);
N = N - 1;
psidot_num = psidot_num*DegElevMatrix(N,N+20);
psidot_den = psidot_den*DegElevMatrix(N,N+20);
psidot_coeff = psidot_num./psidot_den;
psidot_weights = psidot_den; 
psidot = BernsteinPoly(psidot_num,t)./BernsteinPoly(psidot_den,t);


% Plot
figure(4)
plot(t,psidot,'LineWidth',3); hold on
[~,M] = size(psidot_coeff);
plot(linspace(0,tf,M),psidot_coeff,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on


% Subdivide at min 
[minv,iminv] = min(psidot_coeff);
lambdamin = iminv/length(psidot_coeff);
[Cpnum,Pout] = deCasteljau(psidot_num,lambdamin);
[Cpden,Pout] = deCasteljau(psidot_den,lambdamin);
psidot_numA = Cpnum(:, 1:length(psidot_num(1,:)));
psidot_numB = Cpnum(:, length(psidot_num(1,:)):end);
psidot_denA = Cpden(:, 1:length(psidot_num(1,:)));
psidot_denB = Cpden(:, length(psidot_num(1,:)):end);
psidot_coeffA = psidot_numA./psidot_denA;
psidot_weightsA = psidot_denA; 
psidot_coeffB = psidot_numB./psidot_denB;
psidot_weightsB = psidot_denB; 
tA = linspace(0,lambdamin*tf,1000);
tB = linspace(lambdamin*tf,tf,1000);
psidotA = BernsteinPoly(psidot_numA,tA)./BernsteinPoly(psidot_denA,tA);
psidotB = BernsteinPoly(psidot_numB,tB)./BernsteinPoly(psidot_denB,tB);


% Plot
figure(5)
plot(tA,psidotA,'LineWidth',3); hold on
plot(tB,psidotB,'LineWidth',3); hold on
[~,M] = size(psidot_coeffA);
plot(linspace(tA(1),tA(end),M),psidot_coeffA,'o','LineWidth',3);
plot(linspace(tB(1),tB(end),M),psidot_coeffB,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on









return

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D Motion Planning tools
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3D Bernstein Polynomials
p1_N = [1,2,4,5,7; 1,0,2,3,3; 1,1.2,2,3.5,3.5];
p2_N = [1,2,3,5,7; 7,7,6,4,4.5; 1,2,2.5,3,3];
p1N = BernsteinPoly(p1_N,t);
p2N = BernsteinPoly(p2_N,t);

% Plot
figure(1)
plot3(p1N(1,:),p1N(2,:),p1N(3,:),'LineWidth',3); hold on
plot3(p2N(1,:),p2N(2,:),p2N(3,:),'LineWidth',3);
plot3(p1_N(1,:),p1_N(2,:),p1_N(3,:),'o','LineWidth',3);
plot3(p2_N(1,:),p2_N(2,:),p2_N(3,:),'o','LineWidth',3);
set(gca,'fontsize', 40);
axis equal
grid on



%% Speed square
[~,N] = size(p1_N);
N = N - 1;
[tnodes,w,Dm] = LGL_PS(N,tf);
p1_Ndot = p1_N*Dm;
speedsqr_2N = BernsteinProduct(p1_Ndot(1,:),p1_Ndot(1,:))+BernsteinProduct(p1_Ndot(2,:),p1_Ndot(2,:))+BernsteinProduct(p1_Ndot(3,:),p1_Ndot(3,:));
speedsqr2N = BernsteinPoly(speedsqr_2N,t);
speed2N = sqrt(speedsqr2N);


% Plot
figure(2)
plot(t,speed2N,'LineWidth',3); hold on
plot(t,speedsqr2N,'LineWidth',3);
plot(linspace(0,tf,2*N+1),speedsqr_2N,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on




%% Psidot
p1_Nddot = p1_Ndot*Dm;
psidot_num = (BernsteinProduct(p1_Nddot(1,:),p1_Ndot(2,:))-BernsteinProduct(p1_Nddot(2,:),p1_Ndot(1,:)));
psidot_den = speedsqr_2N;
psidot_coeff = psidot_num./psidot_den;
psidot_weights = psidot_den; 
psidot = BernsteinPoly(psidot_num,t)./BernsteinPoly(psidot_den,t);


% Plot
figure(3)
plot(t,psidot,'LineWidth',3); hold on
plot(linspace(0,tf,2*N+1),psidot_coeff,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on


% Degree elevate
[~,N] = size(psidot_num);
N = N - 1;
psidot_num = psidot_num*DegElevMatrix(N,N+10);
psidot_den = psidot_den*DegElevMatrix(N,N+10);
psidot_coeff = psidot_num./psidot_den;
psidot_weights = psidot_den; 
psidot = BernsteinPoly(psidot_num,t)./BernsteinPoly(psidot_den,t);


% Plot
figure(4)
plot(t,psidot,'LineWidth',3); hold on
[~,M] = size(psidot_coeff);
plot(linspace(0,tf,M),psidot_coeff,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on


% Subdivide at min 
[minv,iminv] = min(psidot_coeff);
lambdamin = iminv/length(psidot_coeff);
[Cpnum,Pout] = deCasteljau(psidot_num,lambdamin);
[Cpden,Pout] = deCasteljau(psidot_den,lambdamin);
psidot_numA = Cpnum(:, 1:length(psidot_num(1,:)));
psidot_numB = Cpnum(:, length(psidot_num(1,:)):end);
psidot_denA = Cpden(:, 1:length(psidot_num(1,:)));
psidot_denB = Cpden(:, length(psidot_num(1,:)):end);
psidot_coeffA = psidot_numA./psidot_denA;
psidot_weightsA = psidot_denA; 
psidot_coeffB = psidot_numB./psidot_denB;
psidot_weightsB = psidot_denB; 
tA = linspace(0,lambdamin*tf,1000);
tB = linspace(lambdamin*tf,tf,1000);
psidotA = BernsteinPoly(psidot_numA,tA)./BernsteinPoly(psidot_denA,tA);
psidotB = BernsteinPoly(psidot_numB,tB)./BernsteinPoly(psidot_denB,tB);


% Plot
figure(5)
plot(tA,psidotA,'LineWidth',3); hold on
plot(tB,psidotB,'LineWidth',3); hold on
[~,M] = size(psidot_coeffA);
plot(linspace(tA(1),tA(end),M),psidot_coeffA,'o','LineWidth',3);
plot(linspace(tB(1),tB(end),M),psidot_coeffB,'o','LineWidth',3);
set(gca,'fontsize', 40);
grid on










