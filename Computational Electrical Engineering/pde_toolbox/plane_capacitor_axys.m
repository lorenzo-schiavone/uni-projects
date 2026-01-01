clear; clc; close all

% 6.3 Electrostatic DC in 2D - Plane Capacitor with symmetries

model = createpde();

% Geometry
typ_c=1;x_c=0;y_c=0; R_c=.12; alpha = 2.5; R = R_c*alpha;
gd_c= [typ_c;x_c;y_c;1];
w = 0.2; h=0.01; d=0.02;
typ_1=2; N_1=4; x_1 = w/2*[-1;1;1;-1];y_1 = d/2 + h*[0;0;1;1];
gd_1 = [typ_1;N_1;x_1;y_1]; % first polygon
typ_2=2;N_2=4; x_2 = w/2*[-1;1;1;-1];y_2 = -d/2 -h + h*[0;0;1;1];
gd_2 = [typ_2;N_2;x_2;y_2]; % second polygon
typ_3=2; N_3=4; x_3 = w/2*[-1;1;1;-1]; y_3 = d/2* [-1;-1;1;1];
gd_3 = [typ_3; N_3;x_3;y_3]; % dielectric

gd = zeros(length(gd_1),5); 
gd(1:length(gd_c),1)=gd_c;
gd(1:length(gd_1),2)=gd_1;
gd(1:length(gd_2),3)=gd_2;
gd(1:length(gd_3),4)=gd_3;

% now exploit symmetries - rect to restrict domain
typ_4=2; N_4=4; x_4 = 10*R*[0;0;1;1]; y_4 = 10*R*[0;1;1;0];
gd_4 = [typ_4;N_4;x_4;y_4];
gd(1:length(gd_3),5)=gd_4;

ns = char('background','plate_up','plate_down','dielectric','rect ')';
sf = '(((background -plate_up)-plate_down)+dielectric)*rect'; % * is the intersection

[g, bt] = decsg(gd,sf,ns); 
geometryFromEdges(model, g);

figure
pdegplot(model,"EdgeLabels","on","FaceLabels","on"); title('geometry')
axis equal

%% MESH 
hmax = 0.01; % mesh size
mesh = generateMesh(model,"Hmax",hmax,'HEdge', {[1,2,3,6,8], 0.0025}, 'GeometricOrder','quadratic');

figure
pdemesh(model);
title('mesh')
axis equal

%% BOUNDARY CONDITIONS 
eps0=8.8541878188e-12; epsr = 4;
specifyCoefficients(model ,"m",0,"d",0,"c",@(location ,state)eps0.*(location.x),"a",0,"f",0,'Face',1);
specifyCoefficients(model ,"m",0,"d",0,"c",@(location ,state)eps0*epsr.*(location.x),"a",0,"f",0,'Face',2);

v0 = 1; v1 = -1;
applyBoundaryCondition(model,"dirichlet", "Edge",[5,6], "u",0);
applyBoundaryCondition(model,"dirichlet", "Edge",[1,2,8], "u",v0);
% no flux is automatically imposed

%% RESULTS
results = solvepde(model); 
u = results.NodalSolution;

figure
pdeplot(model ,"XYData",u)
axis equal
title("Numerical Solution");
xlabel("x")
ylabel("y") 
colormap jet
colorbar;

disp("Capacitance:")
fem_nobc = assembleFEMatrices(model); 
We= 4 * u'*fem_nobc.K*u/2;
C = 2 * We / ((v0-v1)^2); % as W = .5 * C * (dV)^2
disp(C)

