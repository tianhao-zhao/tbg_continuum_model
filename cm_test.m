% cm test 
cutoff = 1.1;
N_B = 25;
w = 0;
theta = 0.1;

cm = CM(theta=theta, cutoff=cutoff, N_B=N_B, w=w);
cm.build_g_network;
fig1 = figure();
ax1 = axes(fig1);
scatter(ax1, cm.g_network1(1, :), cm.g_network1(2, :));
hold (ax1, 'on');
scatter(ax1, cm.g_network2(1, :), cm.g_network2(2, :));
axis (ax1, 'equal');
% disp(cm.H)
cm.plot_3d;
cm.plot_kgammamk;

tic;
Bs = 1:20;
E = cm.calculate_e_b(Bs);
fig_b = figure();
ax_b = axes(fig_b);
fig_zoom_in = figure();
ax_zoom_in = axes(fig_zoom_in);
for i = 1:size(E, 1)
    plot(ax_b, Bs, E(i, :)/cm.electron);
    hold(ax_b, 'on');
end
for i = 1:16
    plot(ax_zoom_in, Bs, E(size(E, 1)/2+i, :)/cm.electron);
    hold(ax_zoom_in, 'on');
end
toc;

% calculate the right bottom of matrix
% extract hintra2
E_intra2 = zeros(2*cm.m2*cm.N_B, numel(Bs));
fig_intra2 = figure();
ax_intra2 = axes(fig_intra2);
for i = 1:numel(Bs)
    cm.build_h_b_intra2(Bs(i));
    H_intra2 = cm.H_B((2*cm.m1*cm.N_B+1):(2*cm.m*cm.N_B),...
        (2*cm.m1*cm.N_B+1):(2*cm.m*cm.N_B));
    E_intra2(:, i) = eig(H_intra2);
end
for i = 1:size(E_intra2, 1)
    plot(ax_intra2, Bs, E_intra2(i, :)/cm.electron);
    hold(ax_intra2, 'on');
end

% calculate only one block of h intra2
E_intra2_block = zeros(2*cm.N_B, numel(Bs));
fig_intra2_block = figure();
ax_intra2_block = axes(fig_intra2_block);
for i = 1:numel(Bs)
    cm.build_h_b_intra2(Bs(i));
    H_intra2_block = cm.H_B((2*cm.m1*cm.N_B+1):(2*(cm.m1+1)*cm.N_B),...
        (2*cm.m1*cm.N_B+1):(2*(cm.m1+1)*cm.N_B));
    E_intra2_block(:, i) = eig(H_intra2_block);
end
for i = 1:size(E_intra2_block, 1)
    plot(ax_intra2_block, Bs, E_intra2_block(i, :)/cm.electron);
    hold(ax_intra2_block, 'on');
end