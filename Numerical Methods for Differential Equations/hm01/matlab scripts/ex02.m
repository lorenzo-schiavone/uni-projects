%clear; close all; clc;
disp("-----")
disp("EX02")
fun = @(t,y) -10 * y.^2;
y_0 = 1;
T = 2;
exact_sol = @(y) 1 ./ (1+ 10 * y);

hs = 2.^(-5:-1:-10);
steps = 2.^(5:1:10);
errors = zeros(1,6);
for i = 1:length(hs)
    y_est = runge_kutta4(fun,hs(i), T, y_0);
    errors(i)=abs(y_est(end)-exact_sol(T));
end

figure(3)
loglog(steps, errors, 'b*-')
title('Error for increasing number of steps')
xlabel('number of steps')
ylabel('error at T=2')
hold off

%% display table 
n_steps= steps';
error = errors';
display(table(n_steps, error))

log_steps = log2(steps);
log_error = log2(errors);
slopes = (log_error(2:end) - log_error(1:end-1)) ./ (log_steps(2:end) - log_steps(1:end-1));

fprintf('Slope: %1.4f\n', mean(slopes));
