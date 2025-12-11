% Plane inductor optimization: parameter optimization of planar inductor to maximize the inductance of the device.
lb = [0.0005; 0.00015; 0.0005];
ub = [0.0015; 0.00035; 0.0015];
options = optimoptions('ga','Display','iter', 'PlotFcn','gaplotbestf','MaxTime', 600); %maxtime: 10 minutes
funcipt = @(x) inductanceCompute(x(1),x(2),x(3));
opt = ga(funcipt, 3, [],[],[],[],lb,ub,[],options);