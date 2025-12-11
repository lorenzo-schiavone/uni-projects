% 8.1 - Direct and Iterative Solver
clc;clear;close all
% FROM EXAMPLE 6.1
load("A.mat");
load("b.mat");
load("x.mat");

norm_b = norm(b);
tol = 1e-6; maxit = 200;
Uchol = ichol(A);
[L,U] = ilu(A);

%pcg
[x1, ~,~,~,resvec1] = pcg(A,b,tol,maxit);
[x2, ~,~,~,resvec2] = pcg(A,b,tol,maxit, Uchol, Uchol');
[x3, ~,~,~,resvec3] = pcg(A,b,tol,maxit, L, U);

% gmres
[x4, ~,~,~,resvec4] = gmres(A,b,[],tol,maxit);
[x5, ~,~,~,resvec5] = gmres(A,b,[],tol,maxit, Uchol, Uchol');
[x6, ~,~,~,resvec6] = gmres(A,b,[],tol,maxit, L, U);

figure
semilogy(0:length(resvec5)-1 , resvec5/norm_b, "*-","DisplayName", "cholGMRES");
hold on
semilogy(0:length(resvec2)-1 , resvec2/norm_b, "*-", "DisplayName", "cholCG");
hold on
semilogy(0:length(resvec1)-1 , resvec1/norm_b, "o-", "DisplayName", "CG");
hold on
semilogy(0:length(resvec6)-1 , resvec6/norm_b, "o-","DisplayName", "iluGMRES");
hold on
semilogy(0:length(resvec3)-1 , resvec3/norm_b, "o-","DisplayName", "iluCG");
hold on
semilogy(0:length(resvec4)-1 , resvec4/norm_b, "o-","DisplayName", "GMRES");
legend

%% STEEPEST DESCENT
function [x, resvec] = steepDesc(A, b, x0, tol, maxit)
    x = x0;
    r = b - A*x;
    resvec = zeros(maxit, 1);
    resvec(1) = norm(r);
    it = 1;
    while resvec(it)>tol
        z = A*r;
        alpha = (r'*r)/(r'*z);
        x = x + alpha*r;
        r = r - alpha*z;
        it = it +1;
        resvec(it) = norm(r);
    end
    resvec = resvec(1:it);
end

[x7, resvec7] = steepDesc(A, b, zeros(size(A,1),1), tol, maxit);

figure
semilogy(0:length(resvec7)-1 , resvec7/norm_b, "-","DisplayName", "SteepGrad");
legend

errs = 100/norm(x)* [norm(x1-x), norm(x2 - x),norm(x3 - x),norm(x4- x),norm(x5 - x), norm(x6 - x), norm(x7 - x)]';
T = table(errs);
disp(T)