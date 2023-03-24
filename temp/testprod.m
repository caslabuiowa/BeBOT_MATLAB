clear all
N = 100;

x1 = rand(1,N+1);
x2 = rand(1,N+1);

tic
y = BernsteinProduct(x1,x2);
toc

Prod = ProdMatrix(N); % This can be computed off-line as long as we know the order of B-poly

tic
xaug = x1'*x2;
xaug = reshape(xaug',[(N+1)^2,1]);
ynew = Prod*xaug;
toc