%% ------------------------------------------------------------------------ %%
%% Polygon calculation formula£¬calculate Gkernel (Fx,Fz,Fxz,Fzz) and gravity value (Gx,Gz,Gxz,Gzz) %%
% Fz, Fx, Fxz, Fzz are kernel function matrices.
% Gz represents the vertical component of gravity, Gx denotes the horizontal component of gravity.
% Gxz and Gzz are gravity gradient tensors.
% 'mod' is the density contrast model.
% 'nodex' and 'nodez' are the grid node coordinate files, respectively.
%% load paraments
clear;clc;
load('model');
%% Mothed1£ºUnoptimized Algorithm
total_computation_time = 0;
tic;
[Fz_1, Gz_1] = cal_Fz_mothed_1(mod, gobs, nodex, nodez); % Computation of the gravitational vertical component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fx_1, Gx_1] = cal_Fx_mothed_1(mod, gobs, nodex, nodez); % Computation of the gravitational horizontal component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fxz_1, Gxz_1] = cal_Fxz_mothed_1(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gxz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fzz_1, Gzz_1] = cal_Fzz_mothed_1(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gzz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;
fprintf('Total computation time for Gx, Gz, Gxz, and Gzz using Method 1: %.6f seconds.\n\n', total_computation_time);
%% Mothed2£ºOnly the calculation formulas for the vertical and horizontal sides were optimized.
total_computation_time = 0;
tic;
[Fz_2, Gz_2] = cal_Fz_mothed_2(mod, gobs, nodex, nodez); % Computation of the gravitational vertical component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fx_2, Gx_2] = cal_Fx_mothed_2(mod, gobs, nodex, nodez); % Computation of the gravitational horizontal component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fxz_2, Gxz_2] = cal_Fxz_mothed_2(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gxz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fzz_2, Gzz_2] = cal_Fzz_mothed_2(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gzz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;
fprintf('Total computation time for Gx, Gz, Gxz, and Gzz using Method 2: %.6f seconds.\n\n', total_computation_time);
%% Mothed3£ºThe edge indexing was optimized based on mothed2, meaning each edge is calculated only once.
total_computation_time = 0;
tic;
[Fz_3, Gz_3] = cal_Fz_mothed_3(mod, gobs, nodex, nodez); % Computation of the gravitational vertical component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fx_3, Gx_3] = cal_Fx_mothed_3(mod, gobs, nodex, nodez); % Computation of the gravitational horizontal component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fxz_3, Gxz_3] = cal_Fxz_mothed_3(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gxz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fzz_3, Gzz_3] = cal_Fzz_mothed_3(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gzz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;
fprintf('Total computation time for Gx, Gz, Gxz, and Gzz using Method 3: %.6f seconds.\n\n', total_computation_time);
%% Mothed4£ºThe node indexing was optimized based on mothed3, meaning that the arctan and ln terms for each node are calculated only once.
total_computation_time = 0;
tic;
[Fz_4, Gz_4] = cal_Fz_mothed_4(mod, gobs, nodex, nodez); % Computation of the gravitational vertical component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fx_4, Gx_4] = cal_Fx_mothed_4(mod, gobs, nodex, nodez); % Computation of the gravitational horizontal component.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fxz_4, Gxz_4] = cal_Fxz_mothed_4(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gxz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;

tic;
[Fzz_4, Gzz_4] = cal_Fzz_mothed_4(mod, gobs, nodex, nodez); % Computation of the gravitational tensor component Gzz.
elapsed_time = toc;
total_computation_time = total_computation_time + elapsed_time;
fprintf('Total computation time for Gx, Gz, Gxz, and Gzz using Method 4: %.6f seconds.\n\n', total_computation_time);
%% Plot the graph %%
% Coordinate transformation for subsequent plotting.
[x_node] = transform_nodes(nodex,mod); 
[z_node] = transform_nodes(nodez,mod);

% Plot the results: (1) gravity components; (2) gravity gradient tensor; (3) density contrast model.
figure('Color','w');

ax1 = subplot(3,1,1);
yyaxis left;plot(gobs(21:181), Gx_4(21:181), 'b', 'LineWidth', 1.5);
ylabel('G_x (mGal)');
axis([0 20 -60 60]);
yyaxis right;plot(gobs(21:181), Gz_4(21:181), 'r', 'LineWidth', 1.5);
ylabel('G_z (mGal)');
axis([0 20 -60 -20]);
set(ax1, 'FontSize', 12, 'LineWidth', 1);
grid on;

ax2 = subplot(3,1,2);yyaxis left;plot(gobs(21:181), Gxz_4(21:181), 'b', 'LineWidth', 1.5);
ylabel('G_{xz} (E)');
axis([0 20 -400 400]);
yyaxis right;plot(gobs(21:181), Gzz_4(21:181), 'r', 'LineWidth', 1.5);
ylabel('G_{zz} (E)');
axis([0 20 -400 400]);
set(ax2, 'FontSize', 12, 'LineWidth', 1);
grid on;

ax3 = subplot(3,1,3);
mod_line = reshape(mod, [], 1);
patch('XData', x_node','YData', z_node','CData', mod_line,'FaceColor', 'flat','EdgeColor', 'none');colormap(jet);
axis tight;
set(ax3,'YDir','reverse','FontSize',12,'LineWidth',1);
xlabel('X (km)');
ylabel('Z (km)');
cb = colorbar;ylabel(cb, 'g/cm^3');

pos1 = ax1.Position;
pos2 = ax2.Position;
pos3 = ax3.Position;
common_left  = pos3(1);
common_width = pos3(3);
ax1.Position(1) = common_left;
ax1.Position(3) = common_width;
ax2.Position(1) = common_left;
ax2.Position(3) = common_width;
