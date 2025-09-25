function y_est = simpson_method(fun, h, T, y_0, y_1)
    %% no handle of different last h
    N = floor(T/h)+1;
    y_est = zeros(1,N);
    y_est(1) = y_0;
    y_est(2)=y_1; 
    for n = 3:N
        y_est(n) = fzero(@(y) -y + y_est(n-2) + h/3 * (fun((n-2)*h, y_est(n-2)) + 4 * fun(h*(n-1), y_est(n-1)) + fun(h*n, y)), y_est(n-1));
    end
end