clc; close all; clear
K=5; % number of mesh to consider 
errors = zeros(K,1);
hh = zeros(K,1); % mesh size
iterJs = zeros(K,1); % iteration for pcg with Jacobi preconditioner for stationary solution
iterICs = zeros(K,1); % "" with iChol ""
disp("Loading exact solution...")
Ref = load("solRef.dat");
% % xRef = Ref(:,1);
% % yRef = Ref(:,2);
% % uRef = Ref(:,3);
interp = scatteredInterpolant(Ref(:,1), Ref(:,2), Ref(:,3));
disp("Done!")

uP1ref = [0.2434390; 0.6046775; 0.7454968; 0.7751273];
uP2ref = [0.0772287; 0.2751718; 0.4328630; 0.4805514];
uP3ref = [0.0183037; 0.0928241; 0.1716526; 0.2008722];

for nn = 0:K-1
    mesh_num = nn;
    input_dir = strcat('input/mesh', num2str(mesh_num), '/mesh', num2str(mesh_num), '.');
    coord_file = strcat(input_dir, 'coord');
    topol_file = strcat(input_dir, 'topol');
    bound_file = strcat(input_dir, 'bound');
    track_file = strcat(input_dir, 'track');
    trace_file = strcat(input_dir, 'trace');
    fprintf('Mesh: %i\n', mesh_num);
    tic
    %% nodes
    nodes = load(coord_file);
    tot_nodes= size(nodes,1);
    
    %% elements
    elements = load(topol_file);
    tot_el = size(elements,1);
    areas = zeros(tot_el, 1);
    
    %% computing areas of each element
    for i = 1:tot_el
        matrix = zeros(3,3);
        for j =1:3
            matrix(j,:)= [1 nodes(elements(i, j), :)];
        end
        areas(i,1)=det(matrix)/2;
    end

    %% assembly
    disp('Assembly...')
    A = sparse(tot_nodes,tot_nodes);
    M = sparse(tot_nodes,tot_nodes);
    for k = 1:tot_el
        el_nodes = elements(k, :);
        area = areas(k);

        xi = nodes(el_nodes(1),1);
        yi = nodes(el_nodes(1),2);

        xj = nodes(el_nodes(2),1);
        yj = nodes(el_nodes(2),2);

        xm = nodes(el_nodes(3),1);
        ym = nodes(el_nodes(3),2);
        
        bi= yj-ym;
        ci= xm-xj;

        bj=ym-yi;
        cj=xi-xm;

        bm=yi-yj;
        cm=xj-xi;
        b= [bi bj bm];
        c=[ci cj cm];
        Aloc = (b'*b + c'*c )/(4*area);

        Mloc = area/12 * (eye(3) + ones(3,3));
        for i = 1:3
            row = el_nodes(i);
            for j=1:3
                col = el_nodes(j);
                A(row,col)=A(row,col) + Aloc(i,j);
                M(row,col)=M(row,col) + Mloc(i,j);
            end
        end
    end
    disp('Done!')
    model_time=toc;
    
    fprintf('time to create the model: %f s\n', model_time);
    
    t = 0;
    dt = 0.02;
    T = 10;
    theta = .5;
    
    bound= load(bound_file);
    num_steps = ceil(T/dt);
    tol =1e-8;
    U = zeros(tot_nodes, num_steps+1);
    
    % set up for imposing boundary condition
    A4update = A(:,bound(:,1));
    nBound = size(bound,1);
    %A(bound(:,1),:) = 0;
    A(bound(:,1), :) = sparse(size(bound,1), size(A,2));
    %A(:,bound(:,1)) = 0;
    A(:, bound(:,1)) = sparse(size(A,1), size(bound,1));
    
    for i = 1: nBound
        A(bound(i,1),bound(i,1)) = 1;
    end 
    % figure(1)
    % spy(A)

    K1 = M/dt + theta*A;
    K2 = K1-A; %M/dt - (1-theta)A;
    droptol = 1e-4;
    L = ichol(K1, struct('type','ict','droptol',droptol));
    maxit=100;
    f=zeros(tot_nodes,1);
    t=dt; % t= 0 is U(:,1) = 0 and it is already done.
    disp('Solving PDEs...')
    tic
    for k=1:num_steps
        if t<T/2
            bd_cond = 2*t/T * bound(:,2); %% in this case bound(:,2) = d2tmax 
        end
    
        f_new = 0-A4update *bd_cond; % the minus because it went to the rhs
        f_new(bound(:,1)) =bd_cond;
    
        rhs = K2*U(:,k) + .5*(f + f_new); 
    
        [U(:,k+1), ~,~,~,~] = pcg(K1,rhs,tol, maxit, L, L'); % default x0 is the zero vector
        t=t+dt;
        f=f_new;
    end
    timePDEs = toc;
    disp('Done!')
    fprintf('time to solve PDEs: %f s\n', timePDEs);

    time_to_plot = num_steps*[.25 .5 .75 1]+1;
    track= load(track_file);
    disp('Value Table: ')
    uP1 = U(track(1), time_to_plot)';
    uP2 = U(track(2), time_to_plot)';
    uP3 = U(track(3), time_to_plot)';
    T=table(uP1,uP2,uP3 ,uP1-uP1ref,uP2-uP2ref,uP3-uP3ref, 'VariableNames',["u(P1)","u(P2)","u(P3)", "err u(P1)","err u(P2)","err u(P3)"]);
    disp(T)
    
    
    %%%% STATIONARY SOLUTION: 
    % u' = 0 -> Mu' + Au = f
    % Au = f
    disp('Stationary solution...')
    maxit=500;
    droptol=1e-2;
    L = ichol(A, struct('type','ict','droptol',droptol));
    [x,~,~,iterIC,resvecIC] =  pcg(A,f,tol, maxit, L, L');
    U_stationary=x;
    iterICs(nn+1)=iterIC;
    tic
    Jacobi = diag(diag(A));
    [xJ,~,~,iterJ,resvecJ] =  pcg(A,f,tol, maxit*2, Jacobi); %% check if this is right
    timeJ= toc;
    disp('Done!')
    iterJs(nn+1)=iterJ;
    figure(5)
    semilogy(resvecIC/norm(f), '-*','DisplayName',strcat("mesh: ", num2str(mesh_num)))  
    hold on

    figure(6)
    semilogy(resvecJ/norm(f), '-*','DisplayName',strcat("mesh: ", num2str(mesh_num)))
    hold on

    %%% errors computation
    disp('Computing total error at final time...')
    
    solintp = zeros(tot_nodes);
    for j =1:tot_nodes
        solintp(j) = interp(nodes(j,:));
    end
    err_sq=0;
    for i = 1:tot_el
        id_node = elements(i,:);
        area = areas(i);
        for j = 1:3
            err_sq = err_sq + area/3*(U_stationary(id_node(j))-solintp(id_node(j)))^2;
        end
    end
    errors(nn+1,1)=sqrt(err_sq);
    disp('Done!')
    
    %%% compute hk= max(distanza fra nodi in un elemento)
    hk = 0;
    for i = 1:tot_el
        id_node = elements(i,:);
        for j = 1:2
            NjX = nodes(id_node(j),:);
            for k=j:3
                NkX = nodes(id_node(k),:);
                dist = sum((NkX-NjX).^2);
                if dist>hk
                    hk=dist;
                end
            end 
        end
    end
    hh(nn+1) = sqrt(hk); 
    disp('-------------')
    save(strcat("result_mesh",num2str(mesh_num),".mat"),"U")
end

rk= zeros(K,1);
for k =1:K-1
    rk(k+1) = errors(k)/errors(k+1)*(hh(k+1)/hh(k))^2;
end
disp("FEM convergence summary table: ")
T2 = table(hh, errors, rk);
disp(T2)

disp("Iteration number PCG stationary solution: ")
T3= table(iterJs, iterICs, 'RowNames', ["0", "1", "2", "3", "4"]);
disp(T3)

fig5 = figure(5);
title('Convergence Plot Stationary Solution - IChol')
legend
saveas(fig5,"ConvIChol.jpg")
fig6= figure(6);
title('Convergence Plot Stationary Solution - Jacobi')
legend
saveas(fig6,"ConvJacobi.jpg")