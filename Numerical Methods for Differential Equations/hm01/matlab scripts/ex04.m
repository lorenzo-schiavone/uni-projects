%clc; clear; close all;
disp("-----")
disp("EX04")
nx = 60;
G= numgrid('S',nx);
A = delsq(G)*(nx-1)^2;

y0 = ones(size(A,1), 1); 
fun = @(t,y) -A*y;
T = 0.1;

% y_exact = expm(-T*A)*y0;
% save("y_exact.mat", "y_exact")
load("y_exact.mat");

%%%%%%%%%%
%%%%%%%%%%Interval of absolute stability for RK4  %%%%%%%%%%
max_lambda = eigs(A,1,'lm'); % max eigenvalues

% axisbox = [-4 2 -4 4];
% xa = axisbox(1); xb = axisbox(2); ya = axisbox(3); yb = axisbox(4);
% nptsx = 101;
% nptsy = 101;
% 
% x = linspace(xa,xb,nptsx);
% y = linspace(ya,yb,nptsy);
% [X,Y] = meshgrid(x,y);
% Z = X +1i*Y;

stab_funct = @(h) abs((1/24)*h.^4 + (1/6)*h.^3 + (1/2)*h.^2 + h + 1) - 1;

% val = stab_funct(Z); has real roots near -3 and 0

root1 = fzero(stab_funct, -3); 
root2 = fzero(stab_funct, -1);
 
xP = [root1, root2];
yP = [0, 0];
%%% PLOT REGION OF ABSOLUTE STABILITY RK4
% figure()
% contour(x,y,val,[0 0],'-');
% grid on; 
% hold on
% 
% plot(xP, yP, 'o', 'MarkerSize', 8,  'MarkerFaceColor', 'black');
% hold on;
% plot(xP, yP, '-');
% hold on;
% 
% text(xP(1) + 0.1, yP(1) + .2, sprintf('%.3f', xP(1)), 'FontSize', 12);
% text(.1,.2,sprintf('%1.f', 0), 'FontSize', 12);
% title("|(1/24)*h^4 + (1/6)*h^3 + (1/2)*h^2 + h + 1| = 1")
% hold off

max_norm_root = max(abs(root1), abs(root2));

h_max = max_norm_root/abs(max_lambda);
fprintf('h for RK4 absolute stability: %.4e\n', h_max)

CPUtimes=[];
err_inf = [];
%err_2 = [];

%%% ODE45
disp('ODE45:')
tic
[tt, y_ode45] = ode45(fun, [0 T], y0);
CPUtimes(end+1) = toc;
y_ode45 = y_ode45';
disp('done')

err = abs(y_ode45(:, end) - y_exact);
err_inf(end+1) = max(err);
%err_2(end+1) = norm(err);

%%%CRANK NICHOLSON
disp('CN:')
maxit=40;
droptol = 1e-2;
hh= 10.^-(2:1:5);
for i = 1:length(hh)
    h = hh(i);
    N = floor(T/h);
    tol = h^3;
    tic
    Bh = speye(size(A,1)) + (h/2) * A;
    Ch = speye(size(A,1)) - (h/2) * A;
    L = ichol(Bh, struct('type','ict','droptol',droptol));
    
    y_estCN = y0;
    for n = 1:N
        b = Ch * y_estCN;
        [x,~,~,iter] = pcg(Bh, b, tol, maxit, L, L');
        y_estCN = x ;
    end

    h_last = T - N*h;
    if h_last>0

        Bh_last = speye(size(A,1))+ (h_last/2) * A;
        Ch_last = speye(size(A,1))- (h_last/2) * A;

        b_last = Ch_last * y_estCN;
        [x,~,~,iter] = pcg(Bh_last, b_last, tol, maxit);
        y_estCN = x;
    end
    CPUtimes(end+1) =toc;
    err = abs(y_estCN - y_exact);
    err_inf(end+1) = max(err);
    %err_2(end+1) = norm(err);
    fprintf('h = %i done\n', h)
end

%%%%%%%%%%%%%%% BDF3 %%%%%%%%%%%%%%%%%%%%%%%
disp('BDF3:')
for i = 1:length(hh)
    h = hh(i);
    N = floor(T/h);
    y_estBDF3 = zeros(size(A,1), N);
    y_estBDF3(:,1) = y0;
    %y_estBDF3(:, 2) = expm(-h*A)*y0; % too expensive operation
    %y_estBDF3(:, 3) = expm(-2*h*A)*y0; %too expensive operation
    tol = h^3;
    tic
    %% STEP 1
    tmp1=runge_kutta4(fun, min(1e-4, h/2), h, y_estBDF3(:, 1)); 
    y_estBDF3(:, 2) = tmp1(:,end);
    %% STEP 2
    tmp2=runge_kutta4(fun, min(1e-4, h/2), h, y_estBDF3(:, 2));
    y_estBDF3(:, 3) = tmp2(:,end);

    Bh = speye(size(A,1)) + 6/11 * h * A;
    L = ichol(Bh, struct('type','ict','droptol',droptol));
    for n = 3: N
        b = (18 * y_estBDF3(:,n) - 9 * y_estBDF3(:, n-1) +2 * y_estBDF3(:, n-2))/11;
        [x, ~,~, iter]=pcg(Bh, b, tol, maxit, L, L');
        y_estBDF3(:, n+1) = x;
    end
    CPUtimes(end+1) =toc;
    fprintf('h = %i done\n', h)
    err = abs(y_estBDF3(:,end) - y_exact);
    err_inf(end+1) = max(err);
    %err_2(end+1) = norm(err);
end

methods=["ode45"  repmat("CN", 1, length(hh)) repmat("BDF3", 1, length(hh))];
steps= [length(tt)  T./hh T./hh];

methods = methods'; steps= steps'; err_inf = err_inf';CPUtimes= CPUtimes';%err_2 = err_2';

T = table(methods, steps,CPUtimes, err_inf); %err_2
disp(T)