alpha = 0.2; %20.0;
beta = 0.01; %1.0;
gamma = 0.004; %1.0;
delta = 0.07; %30.0;

x0 = 19; %8;
y0 = 22; %12;
t0 = 0;
T = 300;
h = 1.0e-3;

fun = @(t, z) [z(1)*(alpha - beta*z(2)); z(2)*(gamma*z(1) - delta)];
y_estLV = runge_kutta4(fun, h, T, [x0;y0]);
% 
% N = floor(T/h);
% y_est = zeros(2, N);
% y_est(:,1) = [x0;y0];
% for n = 1:N
%     t = (n-1)* h;
%     k1 = fun(t, y_est(:, n));
%     k2 = fun(t + h/2, y_est(:, n) +1/2 * h * k1);
%     k3 = fun(t + h/2, y_est(:,n) +1/2 * h * k2);
%     k4 = fun(t +h, y_est(:,n) + h*k3);
%     y_est(:, n+1) = y_est(:,n) + 1/6 * h *(k1+ 2*k2 +2*k3 +k4);
% end
% 
% h_last = T - N*h;
% if h_last>0
%     t=N*h;
%     k1 = fun(t, y_est(:, end));
%     k2 = fun(t + h_last/2, y_est(:, end) +1/2 * h_last * k1);
%     k3 = fun(t + h_last/2, y_est(:,end) +1/2 * h_last * k2);
%     k4 = fun(t +h_last, y_est(:,end) + h_last*k3);
%     y_est = [y_est y_est(:,end) + 1/6 * h_last *(k1+ 2*k2 +2*k3 +k4)];
% end

figure(4)
plot(y_estLV(1, :))
hold on
plot(y_estLV(2, :))
hold on
legend('prey', 'predator')
title('Lokta Volterra: prey-predator model')