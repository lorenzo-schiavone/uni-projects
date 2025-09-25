function [x, resvec, iter] = mypcg(A, b, tol, maxit, L)
x = zeros(size(b));
r = b - A *x;
resvec = zeros(maxit);
resvec(1)=norm(r);
p = L' \ (L \ r);
g=p;
rho = r'*g;
k=0;
exit_test= tol*norm(b);
while ((resvec(k+1)>exit_test) && (k<maxit))
    z = A*p;
    alpha = rho/(z'*p);
    x = x + alpha*p;
    r = r -alpha*z;
    g= L' \ (L \ r);
    rho_new = r'*g;
    beta = rho_new/rho;
    p = g + beta*p;
    k=k+1;
    rho = rho_new;
    resvec(k+1)=norm(r);
end
resvec = resvec(1:(k+1));
iter = k;