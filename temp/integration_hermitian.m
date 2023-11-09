%% This function performs hermitian interpolation using Bernstein integration

clear all 
close all
index = 0;
for k = 4:2:200

pause(0.35)

index = index + 1;

tic 

tf = 2.5;
x = 10*rand(k,1)-5;
[Poly, Transf] = HermitianInterpolation(x,tf);

Dm = BernsteinDifferentiationMatrix(k-1,1)';
[n,m] = size(Dm);
midIdx = ceil(m/2);
Dm(midIdx, :) = []; 
Transform = [1 zeros(1,m-1); Dm ; zeros(1,m-1) 1];
CondTransf(index) = cond(Transform);

Dm = BernsteinDifferentiationMatrix(k-1,1)';
[n,m] = size(Dm); 
Transform2 = [1 zeros(1,m-1); Dm];
CondTransf2(index) = cond(Transform2);





figure(1)
subplot(1,2,1)
time = linspace(0,tf,1000);
plot(linspace(0,tf,length(Poly)),Poly,'o','color','b')
hold on
plot(time,BernsteinPoly(Poly',time))
hold off






x = x';
N = length(x)-1;

[~,~,Dm] = BeBOT(N,tf);

BN0 = bernsteinMatrix(N,0);
BN1 = bernsteinMatrix(N,1);

BC = zeros(N+1);
BC(1,:) = BN0;
BC(end,:) = BN1;

j = 1;
for i=2:(N+1)/2
    BC(i,:) = BN0*Dm'^j;
    j = j+1;
end
j = (N+1)/2-1;
for i=(N+1)/2+1:N
    BC(i,:) = BN1*Dm'^j;
    j = j-1;
end

T = inv(BC);
Cp = T*x';

toc

CondT(index) = cond(T);

figure(1)
subplot(1,2,2)
time = linspace(0,tf,1000);
plot(linspace(0,tf,length(Cp)),Cp,'o','color','r')
hold on
plot(time,BernsteinPoly(Cp',time))
hold off

k

end

figure(3)
semilogy([4:2:200],CondT);
figure(4)
semilogy([4:2:200],CondTransf);

