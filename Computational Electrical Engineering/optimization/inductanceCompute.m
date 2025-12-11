% 6.6 Planar Inductor magnetostatic 2d-axisymmetric
function L = inductanceCompute(a,b, d)

model = createpde();
%% GEOMETRY
w=0.015; h = 0.001;
mur = 140; J_source = 1/(a*b);

typ=1; x=0;y=0;bigR = 0.05;
gd_big = [typ;0;0;bigR];

typ=2; N=4; x_1 = [0;w;w;0]; y_1 = [0;0;h;h];
gd_1 = [typ;N;x_1;y_1];

x_2 = [d;d+a;d+a;d]; y_2 = [h;h;h+b;h+b];
gd_2 = [typ;N;x_2;y_2];

typ_2=2; N=4; x_1 = bigR*[0;1;1;0]; y_1 = bigR*[-1;-1;1;1];
gd_4 = [typ_2;N;x_1;y_1];

gd = zeros(length(gd_1),4); 
gd(1:length(gd_big),1)=gd_big;
gd(1:length(gd_1),2)=gd_1;
gd(1:length(gd_2),3)=gd_2;
gd(1:length(gd_4),4)=gd_4;

ns = char('background','core','coil','rect')';
sf = '((core+coil)+background)*rect'; 

[g, bt] = decsg(gd,sf,ns); 
geometryFromEdges(model, g);

%% MESH
hmax = 0.002; % mesh size
mesh = generateMesh(model,"Hmax",hmax, 'Hface', {3, 0.00005, 2, 0.0001},'GeometricOrder','linear');


%% BOUNDARY CONDITION
seps = 1e-8;
mu0=4*pi*1e-7;
specifyCoefficients(model ,"m",0,"d",0,"c",@(location, state) 1./(mu0.*(location.x + seps)) ,"a",0,"f",0,'Face',1);
specifyCoefficients(model ,"m",0,"d",0,"c",@(location, state) 1./(mu0.*(location.x + seps)) ,"a",0,"f",J_source,'Face',3);
specifyCoefficients(model ,"m",0,"d",0,"c",@(location, state) 1./((mur*mu0).*(location.x + seps)),"a",0,"f",0,'Face',2); 

applyBoundaryCondition(model,"dirichlet", "Edge",[9,10,11,12,13], "u",0);

%% SOLUTION AND PLOTS
results = solvepde(model); 
u = results.NodalSolution;

Aphi=u./(model.Mesh.Nodes(1,:)'+seps);

%% INDUCTANCE
t=model.Mesh.Elements; t(end+1,:)=0;

pp=[d+a/2,h+b/2];
F = pdeInterpolant(model.Mesh.Nodes,t,Aphi);

L = - 2*pi*pp(1)*evaluate(F,pp(1),pp(2)); % - so that ga minimize it for getting the max inductance
end
