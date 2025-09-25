clear; clc; close all;
mu=1; dt=1; cmass = 1/dt;

NN = 4;
errs_u1 = zeros(NN,1);
errs_u2 = zeros(NN,1);
errs_p  = zeros(NN,1);
hh = zeros(NN,1);

restart=50; tol=1e-9; maxit=20; 
% options.droptol = 1e-3;
% options.type = "crout";
options.udiag = 1;
options.type='nofill';
options.milu = 'off';

% exact solution
u1 = @(x,y) -cos(x*2*pi).*sin(y*2*pi)+sin(y*2*pi);
u2 = @(x,y) sin(2*pi*x).*cos(2*pi*y)-sin(2*pi*x);
p = @(x,y) 2*pi*(cos(2*pi*y)-cos(2*pi*x));
f1 = @(x,y) cmass * u1(x,y)-4*mu*pi^2*sin(2*pi*y).*(2*cos(2*pi*x)-1)+4*pi^2*sin(2*pi*x);
f2 = @(x,y) cmass * u2(x,y)+4*mu*pi^2*sin(2*pi*x).*(2*cos(2*pi*y)-1)-4*pi^2*sin(2*pi*y);


parfor meshnum = 1:NN
    
    bdnodes_2h = load("mesh/mesh"+num2str(meshnum-1)+"/dirtot.dat");
    nodes_2h = load("mesh/mesh"+num2str(meshnum-1)+"/xy.dat");
    elements_2h = load("mesh/mesh"+num2str(meshnum-1)+"/mesh.dat");
    n_nodes_2h = size(nodes_2h, 1);
    n_el_2h = size(elements_2h,1);
    child = zeros(n_el_2h,4); % map element -> child element index --- maybe it can be not used with idx_el*4+1:4 
    
    n_el_h = 4*n_el_2h;
    elements_h = zeros(n_el_h, 4); 
    nodes_h = nodes_2h(:,:); % to enlarge
    curr_node_idx = n_nodes_2h + 1;
    curr_element_idx = 1;
    
    % map: sorted edge -> index of midpoint. empty if not already done
    midpointEdgeIdx = sparse(n_nodes_2h,n_nodes_2h);
    
    for i = 1:n_el_2h
        el = elements_2h(i,1:3);
        % nodes to add 
        n1 = nodes_2h(el(1),:); n2 = nodes_2h(el(2),:); n3 = nodes_2h(el(3),:);
        new1 = (n1 + n2)/2; new2 = (n2 + n3)/2; new3 = (n1+n3)/2;
        % check if already there, if not add to nodes_h. else pick the index
        if midpointEdgeIdx(el(1), el(2))~= 0
            idx_new1 = midpointEdgeIdx(el(1), el(2));
        else
            nodes_h(curr_node_idx,:) = new1; idx_new1 = curr_node_idx; curr_node_idx = curr_node_idx+1;
            midpointEdgeIdx(el(1), el(2)) = idx_new1; midpointEdgeIdx(el(2), el(1)) = idx_new1; 
        end
        if midpointEdgeIdx(el(2), el(3))~= 0
            idx_new2 = midpointEdgeIdx(el(2), el(3));
        else
            nodes_h(curr_node_idx,:) = new2; idx_new2 = curr_node_idx; curr_node_idx = curr_node_idx+1;
            midpointEdgeIdx(el(2), el(3)) = idx_new2; midpointEdgeIdx(el(3), el(2)) = idx_new2;
        end
        if midpointEdgeIdx(el(3), el(1))~= 0
            idx_new3 = midpointEdgeIdx(el(3), el(1));
        else
            nodes_h(curr_node_idx,:) = new3; idx_new3 = curr_node_idx; curr_node_idx = curr_node_idx+1;
            midpointEdgeIdx(el(3), el(1)) = idx_new3; midpointEdgeIdx(el(1), el(3)) = idx_new3;
        end
        % add elements   
        elements_h(curr_element_idx,:) = [el(1) idx_new1 idx_new3 1];    child(i,1) = curr_element_idx; curr_element_idx = curr_element_idx+1;
        elements_h(curr_element_idx,:) = [el(2) idx_new2 idx_new1 2];    child(i,2) = curr_element_idx; curr_element_idx = curr_element_idx+1;
        elements_h(curr_element_idx,:) = [el(3) idx_new3 idx_new2 3];    child(i,3) = curr_element_idx; curr_element_idx = curr_element_idx+1;
        elements_h(curr_element_idx,:) = [idx_new1 idx_new2 idx_new3 0]; child(i,4) = curr_element_idx; curr_element_idx = curr_element_idx+1;
    end
    n_nodes_h = size(nodes_h,1);

    % computing areas of each element
    areas_2h = zeros(n_el_2h, 1);
    for i = 1:n_el_2h
        matrix = zeros(3,3);
        for j =1:3
            matrix(j,:)= [1 nodes_2h(elements_2h(i, j), :)];
        end
        areas_2h(i,1)=det(matrix)/2;
    end
    
    areas_h = zeros(n_el_h, 1);
    for i = 1:n_el_h
        matrix = zeros(3,3);
        for j =1:3
            matrix(j,:)= [1 nodes_h(elements_h(i, j), :)];
        end
        areas_h(i,1)=det(matrix)/2; 
    end
    
    disp('Assembly...')
    A = sparse(n_nodes_h,n_nodes_h);
    M = sparse(n_nodes_h,n_nodes_h);
    B1 = sparse(n_nodes_2h,n_nodes_h);
    B2 = sparse(n_nodes_2h,n_nodes_h);
    
    for k = 1:n_el_2h % parent element
        el_nod_2h = elements_2h(k, 1:3);
        area_2h = areas_2h(k);
    
        for ii = 1:4 % child element
            el_nod_h = elements_h(child(k,ii),1:3);
            parent_idx = elements_h(child(k,ii),4);

            area_h = area_2h/4; % area_h = areas_h(child(k,ii));
    
            xi = nodes_h(el_nod_h(1),1); xj = nodes_h(el_nod_h(2),1); xm = nodes_h(el_nod_h(3),1);
            yi = nodes_h(el_nod_h(1),2); yj = nodes_h(el_nod_h(2),2); ym = nodes_h(el_nod_h(3),2);
            
            % mesh size
            diam = max([(xi-xj)^2 + (yi-yj)^2, (xi-xm)^2 + (yi-ym)^2, (xm-xj)^2 + (ym-yj)^2]);
            if diam>hh(meshnum)
                hh(meshnum)=sqrt(diam);
            end
        
            bi= yj-ym; bj=ym-yi; bm=yi-yj;
            ci= xm-xj; cj=xi-xm; cm=xj-xi;
            
            b= [bi bj bm];
            c= [ci cj cm];

            Aloc = (b'*b + c'*c )/(4*area_h);
            Mloc = area_h/12 * (eye(3) + ones(3,3));
    
            for i = 1:3
                row = el_nod_h(i);
                for j=1:3
                    col = el_nod_h(j);
                    A(row,col) = A(row,col) + Aloc(i,j);
                    M(row,col) = M(row,col) + Mloc(i,j);
                end
            end
    
            % now the Bs
            if parent_idx==0
                scale = [2 2 2];
            else
                scale = [1 1 1];
                scale(parent_idx) = scale(parent_idx)+3;            
            end
            B1loc = 1 / 12 * diag(scale) * repmat(b,3,1);
            B2loc = 1 / 12 * diag(scale) * repmat(c,3,1);
    
            % assembly the Bs
            for i = 1:3
                row = el_nod_2h(i);
                for j=1:3
                    col = el_nod_h(j);
                    B1(row,col) = B1(row,col) + B1loc(i,j);
                    B2(row,col) = B2(row,col) + B2loc(i,j);
                end
            end
        end
    end
    
    disp('Done!')
    K = mu*A+cmass*M;

    Acal = [K,              zeros(size(A)), -B1'     ;
            zeros(size(A)), K,              -B2'     ;
            -B1,             -B2,            zeros(n_nodes_2h,n_nodes_2h)];
    
    % slow: find boundary nodes
    bdnodes_h = [];
    gamma1 = [];
    for i = 1:n_nodes_h
        ni = nodes_h(i,:);
        if abs(ni(2)-1)<1e-15 % y==1
            gamma1(end+1)=i;
            bdnodes_h(end+1)=i;
        elseif abs(ni(1))<1e-15 % x==0
            bdnodes_h(end+1)=i;
        elseif abs(ni(1)-1)<1e-15 % x==1
            bdnodes_h(end+1)=i;
        elseif abs(ni(2))<1e-15 % y==0
            bdnodes_h(end+1)=i;
        end
    end
    totbdnodes_h = length(bdnodes_h);
    totbdnodes_2h = length(bdnodes_2h);

    % set up boundary values
    uy = zeros(totbdnodes_h,1);
    ux = zeros(totbdnodes_h,1);
    vy= 0; % if we want to add y component to bc in Gamma1
    for i = 1:totbdnodes_h
        if any(gamma1(:) == bdnodes_h(i))
            ux(i)=1.;
            uy(i)=vy;
        end
    end
    bdval = [ux;uy;0];
    
    % impose zero pressure on (0,0)
    % find index of 0,0
    pressure_node = 0;
    for i=1:n_nodes_2h
        if norm(nodes_2h(i,:),1) < 1e-10
            pressure_node = i;
            break
        end
    end
    bdnodes_h = bdnodes_h'; % column for indexing

    % enforce bc
    A4update = Acal(:,[bdnodes_h; bdnodes_h + n_nodes_h; 2*n_nodes_h+pressure_node]); 
    Acal([bdnodes_h;bdnodes_h + n_nodes_h; 2*n_nodes_h+pressure_node],:) = sparse(2*totbdnodes_h+1, size(Acal,2));
    Acal(:,[bdnodes_h;bdnodes_h + n_nodes_h; 2*n_nodes_h+pressure_node]) = sparse(size(Acal,1),2*totbdnodes_h+1);
    Acal([bdnodes_h;bdnodes_h+n_nodes_h; 2*n_nodes_h+pressure_node],[bdnodes_h;bdnodes_h+n_nodes_h; 2*n_nodes_h+pressure_node]) = speye(2* totbdnodes_h+1);
    rhs = 0 - A4update*bdval; 
    rhs([bdnodes_h; bdnodes_h+n_nodes_h; 2*n_nodes_h+pressure_node]) = bdval;
    
    % solve the system with ilu gmres
    tic
    [L,U]=ilu(Acal,options);
    sol=gmres(Acal, rhs, restart, tol, maxit, L, U);
    gmrestime = toc;
    fprintf("time for gmres 2 %f\n", gmrestime)
    
    u1_h = sol(1:n_nodes_h);
    u2_h = sol(1+n_nodes_h:2*n_nodes_h);
    p_h = sol(2*n_nodes_h+1:end);
    
    % plot solution
    if meshnum==4
        figure(10*meshnum + 1)
        trisurf(elements_h(:,1:3), nodes_h(:,1), nodes_h(:,2), u1_h, 'EdgeColor', 'None');
        title("u_x");
        xlabel("x")
        ylabel("y")

        figure(10*meshnum +2)
        trisurf(elements_h(:,1:3), nodes_h(:,1), nodes_h(:,2), u2_h, 'EdgeColor', 'None');
        title("u_y");
        xlabel("x")
        ylabel("y")

        figure(10*meshnum + 3)
        trisurf(elements_2h(:,1:3), nodes_2h(:,1), nodes_2h(:,2), p_h, 'EdgeColor', 'None');
        title("p");
        xlabel("x")
        ylabel("y")
    end
    
    %%%%%%%% KNOWING THE EXACT SOLUTION - for CONVERGENCE ANALYSIS
    
    % compute forcing term rhs
    f1rhs = zeros(n_nodes_h,1); % f1rhs, f2rhs are in the mesh T_h
    f2rhs = zeros(n_nodes_h,1);
    for i=1:n_el_h
        el_nodes = elements_h(i, 1:3);
        area = areas_h(i);
        x = nodes_h(el_nodes,1); y = nodes_h(el_nodes,2);
        f1_el = f1(x,y);
        f2_el = f2(x,y);
        
        f1rhs(el_nodes) = f1rhs(el_nodes) + f1_el*area/3; 
        f2rhs(el_nodes) = f2rhs(el_nodes) + f2_el*area/3;
    end
    f = [f1rhs; 
         f2rhs;
         zeros(n_nodes_2h,1)];

    Acal = [K,              zeros(size(A)), -B1';
            zeros(size(A)), K,              -B2';
            -B1 ,          -B2,            zeros(n_nodes_2h,n_nodes_2h)];
    
    % bd condition
    bdval = [zeros(totbdnodes_h,1);
             zeros(totbdnodes_h,1);
             p(nodes_2h(bdnodes_2h,1),nodes_2h(bdnodes_2h,2))];
    
    A4update = Acal(:,[bdnodes_h; bdnodes_h + n_nodes_h; 2*n_nodes_h+bdnodes_2h]); 
    Acal([bdnodes_h; bdnodes_h + n_nodes_h; 2*n_nodes_h+bdnodes_2h],:) = sparse(2*totbdnodes_h + totbdnodes_2h, size(Acal,2));
    Acal(:,[bdnodes_h; bdnodes_h + n_nodes_h; 2*n_nodes_h+bdnodes_2h]) = sparse(size(Acal,1),2*totbdnodes_h + totbdnodes_2h);
    Acal([bdnodes_h; bdnodes_h+n_nodes_h; 2*n_nodes_h+bdnodes_2h],[bdnodes_h; bdnodes_h+n_nodes_h; 2*n_nodes_h+bdnodes_2h]) = speye(2*totbdnodes_h + totbdnodes_2h);
    rhs = f - A4update*bdval; 
    rhs([bdnodes_h;bdnodes_h + n_nodes_h; 2*n_nodes_h+bdnodes_2h]) = bdval; 
    
    % gmres solver
    tic
    [L,U]=ilu(Acal,options); 
    sol2=gmres(Acal,rhs, restart,tol,maxit, L,U);
    gmrestime = toc;
    fprintf("time for gmres 2 %f\n", gmrestime)
    
    
    u1_h = sol2(1:n_nodes_h);
    u2_h = sol2(1+n_nodes_h:2*n_nodes_h);
    p_h = sol2(2*n_nodes_h+1:end);
    
    % plot solution
    if meshnum==4
        figure(meshnum*10 + 4)
        trisurf(elements_h(:,1:3), nodes_h(:,1), nodes_h(:,2),u1_h, 'EdgeColor', 'none');
        title("u_x");
        xlabel("x")
        ylabel("y")

        figure(meshnum*10 + 5)
        trisurf(elements_h(:,1:3), nodes_h(:,1), nodes_h(:,2),u2_h, 'EdgeColor', 'none');
        title("u_y");
        xlabel("x")
        ylabel("y")
        
        figure(meshnum*10 + 6)
        trisurf(elements_2h(:,1:3), nodes_2h(:,1), nodes_2h(:,2), p_h, 'EdgeColor', 'none');
        title("p");
        xlabel("x")
        ylabel("y")
    end
    
    % error computation \|u1h - u1\|_2,\|u2h - u2\|_2,\|ph - p\|_2 
    erru1=0;
    erru2=0;
    errp =0;
    % for velocity T_h
    for i=1:n_el_h
        el_nodes = elements_h(i, 1:3);
        area = areas_h(i);
        x = nodes_h(el_nodes,1); y = nodes_h(el_nodes,2);
        erru1 = erru1 + area/3 * sum((u1(x,y)-u1_h(el_nodes)).^2);
        erru2 = erru2 + area/3 * sum((u2(x,y)-u2_h(el_nodes)).^2);
    end
    % for pressure T_2h
    for i=1:n_el_2h
        el_nodes = elements_2h(i, 1:3);
        x = nodes_2h(el_nodes,1); y = nodes_2h(el_nodes,2);
        errp = errp + areas_2h(i)/3 * sum((p(x,y)-p_h(el_nodes)).^2);
    end
    
    errs_u1(meshnum) = sqrt(erru1);
    errs_u2(meshnum) = sqrt(erru2);
    errs_p(meshnum) =  sqrt(errp);

end

% error lies in a line in log log plot
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