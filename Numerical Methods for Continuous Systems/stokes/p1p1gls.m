%% STOKES PROBLEM
close all; clear; clc;

errs_u1 = zeros(4,1);
errs_u2 = zeros(4,1);
errs_p  = zeros(4,1);
hh = zeros(4,1);

% option for gmres
restart=50; tol=1e-9; maxit=20;
% options.droptol = 1e-2;
% options.type = "crout";
options.udiag = 1;
options.type='nofill';
options.milu = 'off';

% stabilization parameter
gls = true; %false
if gls
    delta = 0.0001;
else
    delta = 0;
end

% parameter
dt = 1;
mu = 1; cmass = 1/dt; 

% display 
plt = true;

% exact solution
u1 = @(x,y) -cos(x*2*pi).*sin(y*2*pi)+sin(y*2*pi);
u2 = @(x,y) sin(2*pi*x).*cos(2*pi*y)-sin(2*pi*x);
p =  @(x,y) 2*pi*(cos(2*pi*y)-cos(2*pi*x));
f1 = @(x,y) cmass * u1(x,y)-4*mu*pi^2*sin(2*pi*y).*(2*cos(2*pi*x)-1)+4*pi^2*sin(2*pi*x);
f2 = @(x,y) cmass * u2(x,y)+4*mu*pi^2*sin(2*pi*x).*(2*cos(2*pi*y)-1)-4*pi^2*sin(2*pi*y);


for meshnum = 1:4
    
    bdnodes = load("mesh/mesh"  +num2str(meshnum) + "/dirtot.dat");
    nodes = load("mesh/mesh"    +num2str(meshnum) + "/xy.dat");
    elements = load("mesh/mesh" +num2str(meshnum) + "/mesh.dat");

    n_nodes= size(nodes,1);
    elements = elements(:,1:3); 
    n_elements = size(elements,1);
    
    areas = zeros(n_elements, 1);
    % computing areas of each element
    for i = 1:n_elements
        matrix = zeros(3,3);
        for j =1:3
            matrix(j,:)= [1 nodes(elements(i, j), :)];
        end
        areas(i,1)=abs(det(matrix)/2);
    end
    
    disp('Assembly...')
    A = sparse(n_nodes,n_nodes);
    M = sparse(n_nodes,n_nodes);
    B1 = sparse(n_nodes,n_nodes);
    B2 = sparse(n_nodes,n_nodes); %only because we are in p1/p1 they have the same size of A

    for k = 1:n_elements
        el_nodes = elements(k, :);
        area = areas(k);
    
        xi = nodes(el_nodes(1),1); xj = nodes(el_nodes(2),1); xm = nodes(el_nodes(3),1);
        yi = nodes(el_nodes(1),2); yj = nodes(el_nodes(2),2); ym = nodes(el_nodes(3),2);
        
        % mesh size
        diam = max([(xi-xj)^2 + (yi-yj)^2, (xi-xm)^2 + (yi-ym)^2, (xm-xj)^2 + (ym-yj)^2]);
        if diam>hh(meshnum)
            hh(meshnum)=sqrt(diam);
        end
    
        bi= yj-ym; bj=ym-yi; bm=yi-yj;
        ci= xm-xj; cj=xi-xm; cm=xj-xi;
        
        b= [bi bj bm];
        c=[ci cj cm];

        Aloc = (b'*b + c'*c )/(4*area);
        Mloc = area/12 * (eye(3) + ones(3,3));
    
        B1loc = repmat(b, 3, 1) / 6;
        B2loc = repmat(c, 3, 1) / 6;
    
        for i = 1:3
            row = el_nodes(i);
            for j=1:3
                col = el_nodes(j);
                A(row,col) = A(row,col) + Aloc(i,j);
                M(row,col) = M(row,col) + Mloc(i,j);
                B1(row,col)= B1(row,col) + B1loc(i,j);
                B2(row,col) = B2(row,col) + B2loc(i,j);
            end
        end
    end
    disp('Done!')

    K = mu*A+cmass*M;
    
    Acal = [K,              zeros(size(A)), -B1'     ;
            zeros(size(A)), K,              -B2'     ;
            -B1,             -B2,           -delta*A];
    

    % BC
    tot_boundnodes = length(bdnodes);
    gamma1 = [];
    for i= 1:tot_boundnodes
        ni = bdnodes(i);
        if abs(nodes(ni,2)-1)<1e-15 % y==1
            gamma1(end+1)=ni;
        end
    end
    
    uy = zeros(tot_boundnodes,1);
    ux = zeros(tot_boundnodes,1);
    vy= 0;
    for i = 1:tot_boundnodes
        if any(gamma1(:) == bdnodes(i))
            ux(i)=1.;
            uy(i)=vy;
        end
    end
    bdval = [ux;uy;0]; % last zero is zero pressure on a node

    % find index of 0,0
    pressure_node = 0;
    for i=1:tot_boundnodes
        if norm(nodes(bdnodes(i),:),1) < 1e-10
            pressure_node = i;
            break
        end
    end
    
    A4update = Acal(:,[bdnodes; bdnodes + n_nodes; 2*n_nodes+pressure_node]); 
    Acal([bdnodes;bdnodes + n_nodes; 2*n_nodes+pressure_node],:) = sparse(2*tot_boundnodes+1, size(Acal,2));
    Acal(:,[bdnodes;bdnodes + n_nodes; 2*n_nodes+pressure_node]) = sparse(size(Acal,1),2*tot_boundnodes+1);
    Acal([bdnodes;bdnodes+n_nodes],[bdnodes;bdnodes+n_nodes]) = speye(2* tot_boundnodes);
    Acal(2*n_nodes+pressure_node,2*n_nodes+pressure_node) = 1; 
    
    rhs = 0 - A4update*bdval;
    rhs([bdnodes; bdnodes+n_nodes; 2*n_nodes+pressure_node]) = bdval;
    
    tic
    % gmres solver
    [L,U]=ilu(Acal,options);
    sol=gmres(Acal,rhs, restart,tol,maxit, L,U);
    gmrestime = toc;
    fprintf("time for gmres %f\n", gmrestime)
    
    u1h = sol(1:n_nodes);
    u2h = sol(n_nodes+1:2*n_nodes);
    ph = sol(2*n_nodes+1:end);

    % solution plot
    if meshnum==4
        figure(10*meshnum + 1)
        trisurf(elements, nodes(:,1), nodes(:,2), u1h, 'EdgeColor', 'None');
        title("u_x");
        figure(10*meshnum +2)
        trisurf(elements, nodes(:,1), nodes(:,2), u2h, 'EdgeColor', 'None');
        title("u_y");
        figure(10*meshnum + 3)
        trisurf(elements, nodes(:,1), nodes(:,2), ph, 'EdgeColor', 'None');
        title("p");
    end
    
    %%%%%%%% KNOWING THE EXACT SOLUTION - for CONVERGENCE ANALYSIS
    f1rhs = zeros(n_nodes,1);
    f2rhs = zeros(n_nodes,1);
    g_stab = zeros(n_nodes,1);
    for i=1:n_elements
        el_nodes = elements(i, :);
        area = areas(i);
        x = nodes(el_nodes,1); y = nodes(el_nodes,2);
        f1_el = 1/3 * sum(f1(x,y));
        f2_el = 1/3 * sum(f2(x,y));

        f1rhs(el_nodes) = f1rhs(el_nodes) + f1_el*area/3; 
        f2rhs(el_nodes) = f2rhs(el_nodes) + f2_el*area/3;

        % if I kept in memory b and c i wouldnt need to recompute
        xi = x(1); xj = x(2); xm = x(3);
        yi = y(1); yj = y(2); ym = y(3);
    
        bi= yj-ym; bj=ym-yi; bm=yi-yj;
        ci= xm-xj; cj=xi-xm; cm=xj-xi;
    
        b= [bi bj bm];
        c= [ci cj cm];
        % up to here - not that many lines
    
        g_stab_loc = [b' c'] * [f1_el; f2_el]/2; 
        % \grad \phi_i = [b_i, c_i]/(2*area) 

        g_stab(el_nodes) = g_stab(el_nodes) + g_stab_loc;
    
    end

    f = [f1rhs; f2rhs; -delta * g_stab];
    bdval = [u1(nodes(bdnodes,1), nodes(bdnodes,2)); u2(nodes(bdnodes,1), nodes(bdnodes,2)); p(nodes(bdnodes,1), nodes(bdnodes,2))];
    Acal = [K,              zeros(size(A)), -B1'     ;
            zeros(size(A)), K,              -B2'     ;
            -B1,             -B2,           -delta*A];
    
    % BC
    A4update = Acal(:,[bdnodes; bdnodes + n_nodes; 2*n_nodes+bdnodes]); 
    Acal([bdnodes;bdnodes + n_nodes; 2*n_nodes+bdnodes],:) = sparse(3*tot_boundnodes, size(Acal,2));
    Acal(:,[bdnodes;bdnodes + n_nodes; 2*n_nodes+bdnodes]) = sparse(size(Acal,1),3*tot_boundnodes);
    Acal([bdnodes;bdnodes+n_nodes; bdnodes + 2*n_nodes],[bdnodes;bdnodes+n_nodes;bdnodes+ 2*n_nodes]) = speye(3* length(bdnodes));
    rhs = f - A4update*bdval;
    rhs([bdnodes;bdnodes + n_nodes; 2*n_nodes+bdnodes]) = bdval; 
    tic
    % gmres solver
    [L,U]=ilu(Acal,options); % preconditioner
    sol2=gmres(Acal,rhs, restart,tol,maxit, L,U);
    gmrestime = toc;
    fprintf("time for gmres 2 %f\n", gmrestime)
    
    % plot solution
    if meshnum==4
        figure(meshnum*10 + 4)
        trisurf(elements, nodes(:,1), nodes(:,2),sol2(1:n_nodes), 'EdgeColor', 'None');
        title("u_x");
        figure(meshnum*10 + 5)
        trisurf(elements, nodes(:,1), nodes(:,2),sol2(n_nodes+1:2*n_nodes), 'EdgeColor', 'None');
        title("u_y");
        figure(meshnum*10 + 6)
        trisurf(elements, nodes(:,1), nodes(:,2), sol2(2*n_nodes+1:end), 'EdgeColor', 'None');
        title("p");
    end
    
    % error computation \|u1h - u1\|_2^2,\|u2h - u2\|_2^2,\|ph - p\|_2^2 
    u1_h = sol2(1:n_nodes);
    u2_h = sol2(1+n_nodes:2*n_nodes);
    p_h = sol2(2*n_nodes+1:end);
    
    erru1=0;
    erru2=0;
    errp =0;
    for i=1:n_elements
        el_nodes = elements(i, :);
        area = areas(i);
        x = nodes(el_nodes,1); y = nodes(el_nodes,2);
        erru1 = erru1 + area/3 * sum((u1(x,y)-u1_h(el_nodes)).^2);
        erru2 = erru2 + area/3 * sum((u2(x,y)-u2_h(el_nodes)).^2);
        errp  = errp  + area/3 * sum((p(x,y)-p_h(el_nodes)).^2);
    end
    
    errs_u1(meshnum) = sqrt(erru1);
    errs_u2(meshnum) = sqrt(erru2);
    errs_p(meshnum)  = sqrt(errp);

end

figure
loglog(hh,errs_u1, "DisplayName", "u1");
hold on
loglog(hh,errs_u2, "DisplayName", "u2");
hold on
loglog(hh,errs_p, "DisplayName", "p");
hold on
legend
xlabel('h')     
ylabel('Error') 
title('Errors vs h (log-log scale)')
pbaspect([1 1 1])
set(gca, 'DataAspectRatio', [1 1 1])
T = table(hh, errs_u1, errs_u2, errs_p);
disp(T)

% estimate convergence rate approx 2
coef = polyfit(log(hh), log(errs_u1), 1);
fprintf("estimated order of convergence for u1: %1.2f\n", coef(1))
coef = polyfit(log(hh), log(errs_u2), 1);
fprintf("estimated order of convergence for u2: %1.2f\n", coef(1))
coef = polyfit(log(hh), log(errs_p), 1);
fprintf("estimated order of convergence for p: %1.2f\n\n", coef(1))

% they approach order 2 asyntotically
ratio_h = (hh(1:end-1)./hh(2:end)).^2;
r_kU1 = errs_u1(2:end)./errs_u1(1:end-1);
r_kU2 = errs_u2(2:end)./errs_u2(1:end-1);
r_kP = errs_p(2:end)./errs_p(1:end-1);

ratio_u1 = r_kU1.*ratio_h;
ratio_u2 = r_kU2.*ratio_h;
ratio_p = r_kP.*ratio_h;

T2 = table(ratio_u1, ratio_u2, ratio_p);
disp(T2)