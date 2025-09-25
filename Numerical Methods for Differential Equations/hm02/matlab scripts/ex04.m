disp("EX04")
A=gallery ('wathen',100,100);

tol = 1e-10;
maxit = 10000;

itersCG = 0;
itersJ = 0;
itersIC0 = 0;
timeCG = 0;
timeJ = 0;
timeIC0 = 0;
numIter = 5; %100
for i=1:numIter
    x= rand(size(A,1),1);
    b= A *x;
    %% CG
    tic
    [~,~,~,iter,resvecCG] = pcg(A, b, tol, maxit);
    time = toc;
    timeCG = timeCG + time;
    itersCG = itersCG + iter;
    
    %% JACOBI
    tic
    D=diag(diag(A));
    [~,~,~,iter,resvecJCG] = pcg(A, b, tol, maxit, D);
    time = toc;
    timeJ = timeJ + time;
    itersJ = itersJ + iter;
    
    %IC0
    tic
    L=ichol(A);
    [~,~,~,iter,resvecIC0] = pcg(A, b, tol, maxit, L, L');
    time = toc;
    timeIC0 = timeIC0 + time;
    itersIC0 = itersIC0 + iter;
end

method = ["CG"; "Jacobi"; "IC0"];
%residues = [resvecCG(end); resvecJCG(end); resvecIC0(end)];
times = [timeCG;timeJ;timeIC0]/numIter;
iters = [itersCG;itersJ;itersIC0]/numIter;
T = table(method, times, iters);
disp(T)

semilogy(resvecCG, '-*');
hold on;
semilogy(resvecJCG, '-*');
hold on;
semilogy(resvecIC0, '-*');
hold on;
legend('CG', 'Jacobi', 'IC(O)')
xlabel('iteration number');
ylabel('residual norm');
title('Ex04')