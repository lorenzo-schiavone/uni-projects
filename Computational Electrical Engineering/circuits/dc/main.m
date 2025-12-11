%% DC Tableau Analysis
close all; clear; clc;

% data of the problem
[nodes,types,labels,vals,~]=readcir('netlist_01.cir');

% incidence matrix
n = max(nodes(:,1)); 
e = size(labels,1);
A = incidenceMatrix(nodes, n, e);

% component matrices
[R, G, s] = componentMatricesDC(types, vals, e);

% assembly
M = [A zeros(n-1,e) zeros(n-1,n-1);
    zeros(e,e) -eye(e) A';
    R G zeros(e,n-1)];

% right hand side and solution
b = [zeros(n-1+e,1); s];
x = M\b;

% display solution
i = full(x(1:e));
v = full(x(e+1:2*e));
p = full(x(2*e+1:end));
disp("Currents: ")
disp(i)
disp("Voltages: ")
disp(v)
disp("Potentials: ")
disp(p)