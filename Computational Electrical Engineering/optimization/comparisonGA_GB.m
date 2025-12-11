% Comparison between Gradient Based optimization and Stochastic
% Optimization for the Rastrigin
clc; clear; close all;

function [f] = fun_Rastrigin(x,n,A)
fvec = x.^2-A*cos(2*pi.*x);
f = sum(fvec) + A*n;
end

function [f,df] = fun_Rastrigin_and_grad(x,n,A)
fvec = x.^2-A*cos(2*pi.*x);
f = sum(fvec) + A*n;
df = 2.*x + 2*A*pi*sin(2*pi.*x);
end

n = 2; A = 1;
funobj  = @(x) fun_Rastrigin(x,n,A);

% lower bound and upper bound column vector for each variable
bound = 5.12;
lb = -bound*[1;1]; ub = -lb;

N = 100; % number of runs for computing statistics
resultGA = zeros(N,1);
options = optimoptions('ga', 'Display', 'off');
for i = 1:N
    xoptGA = ga(funobj, n,[],[],[],[],lb,ub,[],options); 
    resultGA(i) = funobj(xoptGA);
end
statsGA = [mean(resultGA); std(resultGA); min(resultGA); max(resultGA)];

% gradient based
funobjgrad  = @(x) fun_Rastrigin_and_grad(x,n,A);
resultGB = zeros(N,1);
options = optimoptions('fminunc', 'Display', 'off');

for i = 1:N
    x0 = 2*(rand(2,1)-.5)*bound; % starting point
    xoptGB = fminunc(funobjgrad, x0, options);
    resultGB(i) = funobj(xoptGB);
end
statsGB = [mean(resultGB); std(resultGB); min(resultGB); max(resultGB)];

% statistics table
T = table(statsGB, statsGA, 'VariableNames', {'Gradient Based', 'Genetic'}, 'RowNames', {'mean', 'std','min', 'max'});
disp(T)

