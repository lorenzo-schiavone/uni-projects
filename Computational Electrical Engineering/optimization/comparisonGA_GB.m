% Comparison between Gradient Based optimization and Stochastic
% Optimization for the Rastrigin
clc; clear; close all;

function [f] = fun_Rastrigin(x,n,A)
fvec = x.^2-A*cos(2*pi.*x);
f = sum(fvec) + A*n;
% f = 0;
% for i=1:n
%     f = f+x(i)^2-A*cos(2*pi*x(i));
% end
% f = f + A*n; 
end

function [f,df] = fun_Rastrigin_and_grad(x,n,A)
fvec = x.^2-A*cos(2*pi.*x);
f = sum(fvec) + A*n;
df = 2.*x + 2*A*pi*sin(2*pi.*x);
end

n = 2; A = 10;
funobj  = @(x) fun_Rastrigin(x,n,A);

% lower bound and upper bound column vector for each variable
bound = 5.12;
lb = -bound*[1;1]; ub = -lb;

sample_points = 150;
x = linspace(-bound, bound, sample_points);
y = linspace(-bound, bound, sample_points);
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));
for i = 1:length(x)
    for j = 1:length(y)
        Z(i,j) = funobj([X(i,j),Y(i,j)]);
    end
end

figure(1)
surf(X,Y,Z, 'EdgeColor', 'k')
title("Rastrigin function")
xlabel("x")
ylabel("y")
zlabel("f(x,y)")

N = 10; % number of runs for computing statistics

figure(2)
[C, hContour] = contour(X,Y,Z, 5, 'DisplayName', 'Rastrigin function'); 
title("Gradient-Based vs Stochastic Algorithms (" + num2str(N) + " runs)")
hold on


resultGA = zeros(N,1); xoptGA = zeros(N,2);
options = optimoptions('ga', 'Display', 'off');
for i = 1:N
    xoptGA(i,:) = ga(funobj, n,[],[],[],[],lb,ub,[],options); 
    resultGA(i) = funobj(xoptGA(i,:));
end

figure(2)
plot(xoptGA(:,1),xoptGA(:,2),'r*', 'MarkerSize', 6, 'DisplayName', 'Genetic');
hold on
legend

statsGA = [mean(resultGA); std(resultGA); min(resultGA); max(resultGA)];

% gradient based
funobjgrad  = @(x) fun_Rastrigin_and_grad(x,n,A);
xoptGB = zeros(N,2);
resultGB = zeros(N,1);
options = optimoptions('fminunc', 'Display', 'off');

for i = 1:N
    x0 = 2*(rand(2,1)-.5)*bound; % starting point
    xoptGB(i,:) = fminunc(funobjgrad, x0, options);
    resultGB(i) = funobj(xoptGB(i,:));
end

% figure(2)
plot(xoptGB(:,1),xoptGB(:,2),'b*', 'MarkerSize', 4, 'DisplayName', 'Gradient based');
xlim([-bound, bound])
ylim([-bound, bound])
hold off

statsGB = [mean(resultGB); std(resultGB); min(resultGB); max(resultGB)];

% statistics table
disp(">> Statistics:")
T = table(statsGB, statsGA, 'VariableNames', {'Gradient Based', 'Genetic'}, 'RowNames', {'mean', 'std','min', 'max'});
disp(T)

% table of x vector and f(x) for gb and ga: 10 runs
k = 10;
disp(">> Gradient Based:")
Tgb = table(xoptGB(1:k,1),xoptGB(1:k,2), resultGB(1:k), 'VariableNames', {'x_opt', 'y_opt', 'f(x_opt,y_opt)'});
disp(Tgb)

disp(">> Genetic:")
Tga = table(xoptGA(1:k,1),xoptGA(1:k,2), resultGA(1:k), 'VariableNames', {'x_opt', 'y_opt', 'f(x_opt,y_opt)'});
disp(Tga)

% trend of f_best on a stochastic run
options = optimoptions('ga', 'Display', 'off', 'PlotFcn','gaplotbestf');
ga(funobj, n,[],[],[],[],lb,ub,[],options); 