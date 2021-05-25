close all

N = 12;
T = 1;
t = linspace(0,T,N+1);

lam = 40*sin(t);

[a,w,Dm] = BeBOT(N,T);
% [a,w,Dm] = LGL_PS(N,T);

i = 2;

sum1 = 0;
sum2 = 0;

% Dm = -[N/T*eye(N); zeros(1,N)]+[zeros(1,N);N/T*eye(N)]

for k = 1:N
%     sum1 = sum1 + lam(k)*Dm(k,i)*w(i);
%     sum2 = sum2 + lam(k)*Dm(i,k)*w(k);
    sum1 = sum1 + lam(k)*Dm(k,i);
    sum2 = sum2 + lam(k)*Dm(i,k);
end


lamdot1 = lam*Dm;
lamdot2 = lam*(Dm');



figure
plot(t,lamdot1,'Color','b'); hold on
plot(t,-lamdot2)
ylim([-50 50])

elem = N;

for i = 1:N
    if Dm(i,elem) ~= 0
        count = i;
    end
end

Dm(count-2:count,elem)'
Dmt = Dm';
Dmt(count-2:count,elem)'

Dm