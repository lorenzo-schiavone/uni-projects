% 6.7 Thermal problem in 2D - Chip with Heat Sink - DC and Transient
clc;clear;close all;

model = createpde();

%% GEOMETRY
h1 = 0.002; w = 0.01; h2 = 0.01; w2 = 0.002; h3 = 0.005;

typ=2; N=4; x_1 = [0;w;w;0]; y_1 = [0;0;h1;h1]; % base
gd_1 = [typ;N;x_1;y_1];

typ=2; N=4; x_2 = [0;w;w;0]; y_2 = h1 + [0;0;h2;h2]; % heat sink
gd_2 = [typ;N;x_2;y_2];

typ=2; N=4; x_3 = w2 + [0;w2;w2;0]; y_3 = h1+h2-h3 + [0;0;h3;h3]; % rect1 da togliere
gd_3 = [typ;N;x_3;y_3];

typ=2; N=4; x_4 = 3*w2 + [0;w2;w2;0]; y_4 = h1+h2-h3 + [0;0;h3;h3]; % rect2 da togliere
gd_4 = [typ;N;x_4;y_4];

gd = zeros(length(gd_1),3); 
gd(1:length(gd_1),1)=gd_1;
gd(1:length(gd_2),2)=gd_2;
gd(1:length(gd_3),3)=gd_3;
gd(1:length(gd_4),4)=gd_4;

ns = char('base','sink','rect1','rect2')';
sf = '((base+sink)-rect1)-rect2'; 

[g, bt] = decsg(gd,sf,ns); 
geometryFromEdges(model, g);

figure
pdegplot(model,"EdgeLabels","on","FaceLabels","on"); title('geometry')
axis equal

%% MESH
hmax = 0.0005; % mesh size
mesh = generateMesh(model,"Hmax",hmax,'GeometricOrder','linear');

figure
pdemesh(model);
title('mesh')
axis equal

%% BOUNDARY CONDITIONS
z = 0.01; Text = 20; heatconvecitivity = 50; dissPower = 2;
thermalConductivitySilicon = 148; heatCapacitySilicon = 700; heatDensitySilicon = 2330; 
thermalConductivityHeatSink = 200; heatCapacityHeatSink = 900; heatDensityHeatSink = 2400;

% face 2: chip, d c f
% face 1: heat sink
% for steady state d/dt T = 0, so d = 0
specifyCoefficients(model ,"m",0,"d",0,"c", thermalConductivityHeatSink ,"a",0,"f",0,'Face',1);
specifyCoefficients(model ,"m",0,"d",0,"c", thermalConductivitySilicon ,"a",0,"f",dissPower/(w*h1*z),'Face',2);

edges_idx = [ 2, 3, 4, 5, 8, 10, 11, 12, 13, 14, 15]; %   8, 10];

applyBoundaryCondition(model ,"neumann", "Edge", edges_idx ,'q', heatconvecitivity, 'g', heatconvecitivity*Text);

%% RESULTS AND PLOT

fem = assembleFEMatrices(model,'nullspace');
uc=fem.Kc\fem.Fc;
u=fem.B*uc+fem.ud;

figure
pdeplot(model ,"XYData",u) 
axis equal
title("T [°C], Steady State");
xlabel("x")
ylabel("y") 
colormap jet
colorbar;

%% TRANSIENT 
T = 200; f = 1000;
time_window = 0:T;
x0 = u*0;
p = @(t) (1+0.2*(sin(4*pi*f*t).^2)); % term to be multiplied to the constant load vector fem.Fc
%options = odeset('RelTol', 1e-4);
% transient d~=0
specifyCoefficients(model ,"m",0,"d",heatDensityHeatSink * heatCapacityHeatSink,"c", thermalConductivityHeatSink ,"a",0,"f",0,'Face',1);
specifyCoefficients(model ,"m",0,"d",heatDensitySilicon * heatCapacitySilicon,"c", thermalConductivitySilicon ,"a",0,"f",dissPower/(w*h1*z),'Face',2);
fem = assembleFEMatrices(model,'nullspace');

[tt,xx] = ode15s(@(t,x) fem.M \ (-fem.Kc*x+fem.Fc*p(t)), time_window, x0);%, options);
ut = fem.B*xx.' +repmat(fem.ud,1,length(tt));

node_indices = findNodes(model.Mesh, 'nearest', [0.005; 0.005]);
figure
plot(tt,(ut(node_indices,:)))
xlabel('time [s]')
title("Temperature in a point [°C]")
drawnow
