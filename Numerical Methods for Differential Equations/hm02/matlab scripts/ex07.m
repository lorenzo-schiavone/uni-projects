%%% EX07 %%%
disp("EX07")
load mat13041.rig
A = spconvert(mat13041);
x_exact = 1./sqrt(1:size(A,1));
x_exact = x_exact';
b=A*x_exact;
tol= 1e-12;
maxit = 550;
x0=zeros(size(A,1),1);

options.droptol = 1e-2;
options.type = 'crout';
[L,U] = ilu(A, options);
restarts = [10; 20; 30; 50];


N = length(restarts); % number of line in the plot
Legend=cell(N,1);

total_it = [];
res = [];
cputimes = [];
figure(6)
for i = 1:N
    restart = restarts(i,1);
    tic
    [x,flag, relres, iter,resvec ] = gmres(A,b,restart,tol ,maxit,L,U,x0);
    cputime= toc;
    cputimes(end+1) = cputime;
    total_it(end+1) = iter(1)*restart + iter(2);
    res(end+1) = resvec(end);
    semilogy(resvec, '-*')
    hold on
    Legend{i}=strcat("restart: ",num2str(restart));
end
legend(Legend)
xlabel('iteration number');
ylabel('residual norm');
title('Ex07')
cputimes = cputimes';
total_it = total_it';
res = res';

T= table(restarts, total_it, res, cputimes);
disp(T)