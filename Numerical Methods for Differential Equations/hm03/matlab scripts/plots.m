close all
for nn = 0:4
    mesh_num = nn;
    input_dir = strcat('input/mesh', num2str(mesh_num), '/mesh', num2str(mesh_num), '.');
    track_file = strcat(input_dir, 'track');
    trace_file = strcat(input_dir, 'trace');
    U = load( strcat("result_mesh",num2str(mesh_num), ".mat"));
    U=U.U;
    dt = 0.02;
    T = 10;
    num_steps = ceil(T/dt);
    track= load(track_file);
    %%% TRACK PLOTS
    fig=figure(mesh_num+1);
    for id_node = track
        plot(dt*(1:num_steps+1),U(id_node, 1:num_steps+1), 'LineWidth', 2)
        hold on
    end
    xlabel('t')
    ylabel('u')
    ylim([0,1])
    xlim([0,10])
    legend("P1","P2", "P3", 'Location','northwest')
    saveas(fig,strcat("track",num2str(mesh_num), ".jpg"))

    trace= load(trace_file);
    boundary_nodes = trace(:,1);

    time_to_plot = num_steps*[.25 .5 .75 1]+1;
    xx=linspace(0, 6, length(boundary_nodes)); %6 is the length of the external boundary
    fig=figure(mesh_num+11);
    for time = time_to_plot
        plot(xx, U(boundary_nodes, time), 'LineWidth', 2)
        hold on
    end
    ylabel('u')
    xlabel('arc length')
    legend('t = 2.5', 't = 5','t = 7.5','t = 10')
    ylim([0,1])
    saveas(fig,strcat("trace",num2str(mesh_num), ".jpg"))
end

%%% MESH2 TRIANGULATION PLOT
nn=2;
mesh_num = nn;
input_dir = strcat('input/mesh', num2str(mesh_num), '/mesh', num2str(mesh_num), '.');
coord_file = strcat(input_dir, 'coord');
topol_file = strcat(input_dir, 'topol');
track_file = strcat(input_dir, 'track');

track= load(track_file);
nodes = load(coord_file);
elements = load(topol_file);
fig=figure(mesh_num+21);
h= patch('Vertices', nodes, 'Faces', elements, 'FaceVertexCData', zeros(size(nodes,1),1), ...   % Color data at vertices
      'FaceColor', 'interp');
hold on
plot([-1 -.1], [1,1], 'b-', 'LineWidth', 4) %GammaD1
hold on
plot([.1 1], [1,1], 'r-', 'LineWidth', 4) %GammaD2
text(-.5,1.08, "\Gamma_{D1}", 'FontSize',12, 'Color', 'blue')
text(.5,1.08, "\Gamma_{D2}", 'FontSize',12, 'Color', 'red')
hold on
trackp = nodes(track,:);
plot(trackp(:,1),trackp(:,2), '.', 'MarkerSize', 20);
for i=1:3
    text(trackp(i,1),trackp(i,2) +.1, strcat('P', num2str(i)), 'FontSize',12, 'Color', 'yellow')
end

saveas(fig,strcat("mesh",num2str(mesh_num),".jpg"))
%%%%

nn=4;
mesh_num = nn;
input_dir = strcat('input/mesh', num2str(mesh_num), '/mesh', num2str(mesh_num), '.');
coord_file = strcat(input_dir, 'coord');
topol_file = strcat(input_dir, 'topol');
nodes = load(coord_file);
elements = load(topol_file);
U = load( strcat("result_mesh",num2str(mesh_num), ".mat"));
U=U.U;
num_steps = 500;
dt = 0.02;

fig=figure(41);
h= patch('Vertices', nodes, 'Faces', elements, 'FaceVertexCData', U(:,1), 'FaceColor', 'interp', 'EdgeColor', 'none');
colorbar
for k=time_to_plot
    set(h, 'FaceVertexCData', U(:,k))
    title(sprintf('t = %1.1f s', (k-1)*dt));
    clim([0,1])
    colorbar('Limits',[0 1])
    drawnow
    saveas(fig,strcat("sol_mesh",num2str(mesh_num), "_time",num2str((k-1)*dt) ,".jpg"))
end