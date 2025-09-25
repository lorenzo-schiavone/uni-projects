disp("EX03")
n=1e4;
custdiag = [200.*(1:5), ones(1,n - 5)];
A = diag(custdiag);
tol =1e-8;
maxit=50;
b= rand(n,1);

[x, resvec, k] = mypcg(A, b, tol, maxit, eye(n));

figure(2)
semilogy(0:length(resvec)-1,resvec, '*-')
xlabel('iteration number');
ylabel('residual norm');
title('Ex03')