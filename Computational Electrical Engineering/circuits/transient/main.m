%% TRANSIENT
clc; clear; close all;

% data
T = 5; 
[nodes,types,labels,vals,~]=readcir('netlist_05_RL.cir');

% incidence matrix
n = max(nodes(:,1)); 
e = size(labels,1);
A = incidenceMatrix(nodes, n, e);

% component matrices
[R, G, L, C, s] = componentMatrices_nonDC(types, vals, e); 

% assembly
M1 = [A zeros(n-1,e) zeros(n-1,n-1);
    zeros(e,e) -eye(e) A';
    R G zeros(e,n-1)];
M2 = [zeros(n-1+e,2*e+n-1);
      L C zeros(e,n-1)];
b0 = [zeros(n-1+e,1);s];

% IC: inductance -> open switch, i.e. resistor with infty resistance 
idxL = find(types =='L');
M1open = M1(:,:);
M1open(n-1+e+idxL,idxL) = 1e9 * numel(idxL);
x0 = M1open\b0;

% solution with ode23
F = @(t,x) -M1*x + b0;
options = odeset('Mass', M2, 'RelTol', 1e-5);
[time,sol] = ode23t(F,[0,T],x0,options);

figure
plot(time,sol(:,1:e), '-*')
title("Currents")
legend(compose('i%d', 1:e))

figure
plot(time,sol(:,e+1:2*e), '-*')
title("Voltages")
legend(compose('v%d', 1:e))

figure
plot(time,sol(:,2*e+1:end), '-*')
title("Potentials")
legend(compose('p%d', 1:n-1))