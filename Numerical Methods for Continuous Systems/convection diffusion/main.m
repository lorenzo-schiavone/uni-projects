clear; clc; close all;
% mesh
N = 20;    % Number of divisions per edge
h = 1/N * sqrt(2); % visual inspection of the mesh
l = 1/N; % length of edge in the border

elements = load("./mesh/mesh.dat");
nodes = load("./mesh/xy.dat");
dirnodes = load("./mesh/dirnod.dat");
dirval = load("./mesh/dirval.dat");

% parameter
beta = [1;3];
beta_norm = norm(beta);

ii=2;
Ks = [1e-3, 1e-2, 1e-1, 1];
K = Ks(ii);

taus = [0.005 0.01 0 0];
tau = taus(ii); 
stabilized = true;
if ~stabilized
    tau = 0;
end
tot_el = size(elements,1);
tot_nodes = size(nodes, 1);

%% assembly
disp('Assembly...')

A = sparse(tot_nodes, tot_nodes);
B = sparse(tot_nodes, tot_nodes);
S = sparse(tot_nodes, tot_nodes);

tic
for k = 1:tot_el
    matrix = zeros(3,3);
    for j =1:3
        matrix(j,:)= [1 nodes(elements(k, j), :)];
    end
    area=det(matrix)/2; 
    el_nodes = elements(k, :);
    xi = nodes(el_nodes(1),1); xj = nodes(el_nodes(2),1); xm = nodes(el_nodes(3),1);
    yi = nodes(el_nodes(1),2); yj = nodes(el_nodes(2),2); ym = nodes(el_nodes(3),2);

    bi= yj-ym; bj=ym-yi; bm=yi-yj;
    ci= xm-xj; cj=xi-xm; cm=xj-xi;

    b = [bi bj bm]; c = [ci cj cm];
    % stiffness
    Aloc = K * (b'*b + c'*c ) / (4 * area); 
    % convection
    sigma = beta' * [b; c];         
    Bloc = (1 / 6) * (ones(3,1) * sigma); 
    % supg 
    Sloc = tau * h / (4*area*beta_norm*K) * (sigma' * sigma);

    % assembly
    for i = 1:3
        row = elements(k,i);
        for j = 1:3
            col = elements(k,j);
            A(row,col) = A(row,col) + Aloc(i,j);
            B(row,col) = B(row,col) + Bloc(i,j);
            S(row,col) = S(row,col) + Sloc(i,j);
        end
    end
end

disp('Done!')
model_time=toc;

fprintf('time to create the model: %f s\n', model_time);

H = A+B+S;
% enforce bc
H4update = H(:,dirnodes);
nBound = length(dirnodes);
H(dirnodes, :) = sparse(nBound, size(H,2));
H(:, dirnodes) = sparse(size(H,1), nBound);
H(dirnodes(:,1),dirnodes(:,1)) = speye(nBound);
q = 0 - H4update * dirval; 
q(dirnodes) = dirval;

% solve the non symmetric linear system with gmres
restart=10; tol=1e-9; maxit=20;
setup.type='nofill';
setup.milu='off';
[L,U]=ilu(H,setup);
z=gmres(H,q,restart,tol,maxit,L,U);

figure()
trisurf(elements, nodes(:,1), nodes(:,2),z, 'EdgeColor', 'k', 'FaceColor', 'none');
title("K = " + num2str(K) + ", tau = " +num2str(tau))
xlabel("x")
ylabel("y")
zlim([min(z) max(z)])
view(150, 15)
pbaspect([1 1 1])


% mass conservation 
% convective flux is always -beta_x - beta_y*3/10  = -1 - 3*.3 = -1.9
flux_D = A * z;
fprintf("diffusive flux: %.3e\n", sum(flux_D(dirnodes)));
