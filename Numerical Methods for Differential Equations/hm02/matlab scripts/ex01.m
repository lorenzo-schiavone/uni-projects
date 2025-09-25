%%%%%%%% EX01 %%%%%%%%%%
close all; clear; clc;
disp("EX01")
A= delsq(numgrid ('S', 102));

n = size(A,1);
b=A*ones(n,1) ;
tol = 1e-8; maxit = 200;

N = 4; % number of line in the plot
Legend=cell(N,1);
Legend{1} = 'MY CG';
Legend{2} = 'MY IC0';
Legend{3} = 'MATLAB CG';
Legend{4} = 'MATLAB IC(0)';

%myCG
[x1, resvec1, iter1] = mypcg(A, b, tol, maxit, speye(size(A,1)));
%matlab CG
[xM1, flagM1 , relresM1 , iterM1 , resvecM1]=pcg(A, b, tol, maxit) ;

%myPCG
L= ichol(A);
[x2, resvec2, iter2] = mypcg(A, b, tol, maxit, L);

%matlab PCG
%L= ichol(A);
[xM2, flagM2 , relresM2 , iterM2 , resvecM2]=pcg(A, b, tol, maxit, L, L') ;

figure(1)
semilogy(0:iter1 ,resvec1, '*'); hold on; 
semilogy(0:iter2 ,resvec2, '*'); hold on;
semilogy(0:iterM1 ,resvecM1, 'k-', 'LineWidth',1.5); hold on;
semilogy(0:iterM2 ,resvecM2, '-','Color',"#808080" ,'LineWidth',1.5);
xlabel('iteration number'); ylabel ( 'Residual norm' ) ;

legend(Legend, 'FontSize',12);
title('Ex01')


