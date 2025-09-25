clear; close all; clc;
fun = @(t,y) -5 *y;
exact_sol = @(y) exp(-5 .* y);
h = .05;
T = 6;
y0 =1;
xx = 0:h:T;
exact = exact_sol(xx);

%%% y1 as exact soultion
y1_exact = exact_sol(h);

simps_exact = simpson_method(fun, h, T, y0, y1_exact);

%%% y1 with RK4
% t=0;
% k1 = fun(t, y0);
% k2 = fun(t + h/2, y0 +1/2 * h * k1);
% k3 = fun(t + h/2, y0 +1/2 * h * k2);
% k4 = fun(t +h, y0 + h*k3);
% y1_rk4 = y0 + 1/6 * h *(k1+ 2*k2 +2*k3 +k4);
y1_rk4= runge_kutta4(fun,h,h,y0);
simps_rk = simpson_method(fun, h, T, y0, y1_rk4(:,end));

%%% y1 with FE
y1_fe = y0 + h * fun(h, y0);
simps_fe = simpson_method(fun, h, T, y0, y1_fe);


figure(1)
title('Error with different second initial value')
semilogy(xx, abs(exact-simps_exact))
hold on
semilogy(xx, abs(exact-simps_rk))
hold on
semilogy(xx, abs(exact-simps_fe))
hold off 
legend('exact second initial value', 'rk4 second initial value', 'fe second initial value')
xlabel('times')
ylabel('error')

figure(2)
title('Error with different second initial value')
plot(xx, (exact-simps_exact))
hold on
plot(xx, (exact-simps_rk))
hold on
plot(xx, (exact-simps_fe))
hold off 
legend('exact second initial value', 'rk4 second initial value', 'fe second initial value')
xlabel('times')
%ylim([-.2,.2])