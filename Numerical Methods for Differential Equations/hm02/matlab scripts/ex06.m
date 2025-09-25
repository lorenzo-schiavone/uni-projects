disp("EX06")
load mat13041.rig
A = spconvert(mat13041);
x_exact = 1./sqrt(1:size(A,1));
x_exact = x_exact';
b=A*x_exact;
tol= 1e-10;
maxit = 550;
x0=zeros(size(A,1),1);

%% PRECGMRES %%
%disp('myprecgmres')
options.droptol = 1e-1;
options.type = 'crout';
[L,U] = ilu(A, options);

%tic
[xP,iterP ,resvecP ,flagP] = myprecgmres(A,b,tol ,maxit,x0,L,U);
%toc
figure(5)
semilogy(resvecP, '-')
hold on
N = 2; % number of line in the plot
Legend=cell(N,1);
Legend{1}='mypregmres';

% fprintf('error: %.2E\n',norm(abs(xP-x_exact)))
% fprintf('residue: %.2E\n', resvecP(end))
% fprintf('true residue: %.2E\n', norm(b-A*xP))
% fprintf('number of iterations: %5u\n', iterP)
% fprintf('flag: %1u\n', flagP)

%disp('matlab gmres')
%tic
[xM,flagM, relsresM, iterM,resvecM ] = gmres(A,b,maxit,tol ,maxit,L,U,x0);
%toc

semilogy(resvecM, '*')
hold on 
Legend{2}='matlab gmres';

% fprintf('error: %.2E\n',norm(abs(xM-x_exact)))
% fprintf('residue: %.2E\n', resvecM(end))
% fprintf('true residue: %.2E\n', norm(b-A*xM))
% fprintf('number of iterations: %5u\n', iterM(end))
% fprintf('flag: %1u\n', flagM)
xlabel('iteration number');
ylabel('residual norm');
legend(Legend)
title('Ex06')

methods = ["myprecgmres"; "gmres"];
errors = [norm(abs(xP-x_exact));norm(abs(xM-x_exact))];
residues = [resvecP(end);resvecM(end)];
true_res = [norm(b-A*xP);norm(b-A*xM)];
iterations = [ iterP;iterM(end)];

T = table(methods, iterations, errors, residues, true_res);
disp(T)
