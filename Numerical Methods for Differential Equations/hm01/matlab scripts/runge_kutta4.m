function y_est = runge_kutta4(fun, h, T, y_0)
    N = floor(T/h)+1;
    y_est = zeros(length(y_0),N);
    y_est(:,1) = y_0;
    for n = 2:N
        t = n*(h-1);
        k1 = fun(t, y_est(:,n-1));
        k2 = fun(t + h/2, y_est(:,n-1) +1/2 * h * k1);
        k3 = fun(t + h/2, y_est(:,n-1) +1/2 * h * k2);
        k4 = fun(t +h, y_est(:,n-1) + h*k3);
        y_est(:,n) = y_est(:,n-1) + 1/6 * h *(k1+ 2*k2 +2*k3 +k4);
    end
end