disp("EX08")
A= load('ML_laplace.mtx'); A= spconvert(A);
nnzA= nnz(A);
n= size(A,1);
x_exact = ones(n,1);
x0= zeros(n,1);
b=A*x_exact;
tol= 1e-12;
maxit = 550;

droptols = [2e-2; 1e-2; 3e-3; 1e-3; 1e-4; 1e-5];
restart = 50;
options.type = 'crout';

N = length(droptols); % number of line in the plot
Legend=cell(N,1);

t_precs = [];
t_sols = [];
cputimes = [];
total_it = [];
res =[];
rho=[];
figure(7)
for i = 1:N
    options.droptol = droptols(i);
    tic
    [L,U] = ilu(A, options);
    t_precs(end+1)= toc;
    nnzL = nnz(L); nnzU=nnz(U);
    tic
    [x,flag, relres, iter,resvec ] = gmres(A,b,restart,tol ,maxit,L,U,x0);
    t_sols(end+1)= toc;
    cputimes(end+1) = t_sols(end)+ t_precs(end);
    total_it(end+1) = iter(1)*restart + iter(2);
    res(end+1)= resvec(end);
    rho(end+1)= (nnzL +nnzU - n)/nnzA;
    semilogy(resvec, '-*')
    hold on 
    Legend{i} = strcat("droptol: ", num2str(droptols(i)));
end
legend(Legend);
xlabel('iteration number');
ylabel('residual norm');
title('Ex08')

t_precs = t_precs';
t_sols =t_sols';
cputimes = cputimes';
total_it = total_it';
res =res';
rho=rho';

T = table(droptols, total_it,t_precs, t_sols, cputimes, res, rho);
disp(T)