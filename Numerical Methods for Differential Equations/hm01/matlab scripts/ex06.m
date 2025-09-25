disp("-----")
disp("EX06")

fun = @(t,y) -10 * y^2;
y0 = 1;
T =100;
h = 1e-3;

sol= @(t) 1./(1+10*t);

%%%%%%%% ADAPTIVE TS %%%%%%%%%
tol = 1e-8;
times = [0];
y_est = [y0];

t = 0;
disp('Time with adaptive timestep:')
tic
while t< T
    one_stage=fzero(@(y) -y + y_est(end) + h/2* (fun(t, y_est(end)) +fun(t+h, y) ),y_est(end));

    two_stage_half=fzero(@(y) -y + y_est(end) + h/4* (fun(t, y_est(end)) +fun(t+h/2, y) ),y_est(end));
    two_stage = fzero(@(y) -y + two_stage_half + h/4* (fun(t+h/2, two_stage_half) +fun(t+h, y) ),two_stage_half);

    err = abs(two_stage - one_stage);
    t_next = min(t + h * (3/4 * tol/err)^(1/3), T); %% to not overcome T
    h = t_next -t;
    y_est(end+1)= fzero(@(y) -y + y_est(end) + h/2* (fun(t, y_est(end)) +fun(t_next, y) ),y_est(end));
    t = t_next;
    times(end+1)= t;
end
toc
fprintf('number of steps: %5.0f\n', size(times,2)-1);
%%%%%%%% CONSTANT TIME STEP %%%%%%%%%
h = 1e-3;
N = floor(T/h)+1;
y_est_const = zeros(1,N+1);
y_est_const(1) = y0;
fprintf('\n')
disp('Time with fixed timestep:')
tic
for n = 2:N+1
    y_est_const(n) = fzero(@(y) -y + y_est_const(n-1) + h/2* (fun(h*(n-1), y_est_const(n-1)) +fun(h*n, y) ),y_est_const(n-1));
end
toc
fprintf('number of steps: %5.0f\n', size(0:h:N*h,2)-1);

figure(5)
loglog(0:h:N*h, abs(sol(0:h:N*h) - y_est_const), '.')
hold on
loglog(times, abs(sol(times)- y_est), '.')
hold on

title('Errors in absolute value')
legend('fixed timesteps', 'adaptive timesteps')
% ylim([0 2*tol])
% xlim([-1, 100])
xlabel('time')
ylabel('error')


figure(6)
loglog(h*ones(1,N-1),'.')
hold on 
loglog(times(2:end)-times(1:end-1), '.')
hold on 
title('Timestep h at each step')
legend('fixed timesteps', 'adaptive timesteps')
xlabel('step number')
ylabel('h')




