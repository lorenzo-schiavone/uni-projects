disp("EX02")
nx_arr = [102; 202; 402; 802];
droptol_arr = [1e-2 1e-3];
tol = 1e-8;
maxit= 5000;
CGiters = zeros(4,1); %[];
IC0iters = zeros(4,1); %[];
IC1iters = zeros(4,1); %[];
IC2iters = zeros(4,1); %[];
% I recorded also the cputimes for curiosity
CGtimes = zeros(4,1); %[];
IC0times = zeros(4,1); %[];
IC1times = zeros(4,1); %[];
IC2times = zeros(4,1); %[];
for i=1:4
    nx=nx_arr(i);
    A= delsq(numgrid ('S', nx));
    x_exact= 1./sqrt(1:size(A,1));
    b = A* x_exact';
    tic
    [~,~,~,iter,~]= pcg(A, b, tol, maxit);
    time = toc;
    CGiters(i) = iter;
    CGtimes(i) = time;
    
    tic
    L=ichol(A);
    [~,~,~,iter,~] = pcg(A, b, tol, maxit, L, L');
    time = toc;
    IC0iters(i) = iter;
    IC0times(i) = time;

    tic
    L = ichol(A, struct('type','ict','droptol',droptol_arr(1)));
    [~,~,~,iter,~] = pcg(A, b, tol, maxit, L, L');
    time = toc;
    IC1times(i) = time;
    IC1iters(i) = iter;

    tic
    L = ichol(A, struct('type','ict','droptol',droptol_arr(2)));
    [~,~,~,iter,~] = pcg(A, b, tol, maxit, L, L');
    time = toc;
    IC2times(i) = time;
    IC2iters(i) = iter;
end


% nx_arr = nx_arr'; IC2times=IC2times'; IC2iters=IC2iters';IC1times=IC1times';IC1iters= IC1iters';
% IC0iters = IC0iters'; IC0times=IC0times';CGiters=CGiters'; CGtimes= CGtimes';

T=table(nx_arr,CGiters, CGtimes,  IC0iters, IC0times, IC1iters, IC1times, IC2iters, IC2times);
disp(T)

