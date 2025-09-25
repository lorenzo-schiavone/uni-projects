disp("EX05")
load mat13041.rig
A = spconvert(mat13041);
x_exact = 1./sqrt(1:size(A,1));
x_exact = x_exact';
b=A*x_exact;
tol= 1e-10;
maxit = 550;
x0=zeros(size(A,1),1);

%disp('mygmres')
tic
[x, iter, resvec, flag] = mygmres(A,b,tol,maxit,x0);
time = toc;
figure(4)
semilogy(resvec, '*')
hold on
N = 1; % number of line in the plot
Legend=cell(N,1);
Legend{1}='mygmres';
legend(Legend)
xlabel('iteration number');
ylabel('residual norm');
title('Ex05')
fprintf('gmres cpu time: %.3f sec\n', time)
% fprintf('error: %.2E\n',norm(abs(x-x_exact)))
% fprintf('residue: %.2E\n', resvec(end))
% fprintf('number of iterations: %5u\n', iter)
% fprintf('flag: %1u\n', flag)


