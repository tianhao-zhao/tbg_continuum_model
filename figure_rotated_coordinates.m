% linewidth
lw = 1.5;
len = 1;

% K K' points
K = [sqrt(3)/2, 1/2];
Kp = [sqrt(3)/2, -1/2];

fig = figure();
ax = axes(fig);
plot_single_hexagon(ax, [0, 0], 1, pi/6, 0, lw);
add_local_coordinates(ax, K, 0, len, [0.5, 0.5, 0.5], lw, '--');
add_local_coordinates(ax, K, pi/6, len, 'black', lw, '-');
add_local_coordinates(ax, Kp, 0, len, [0.5, 0.5, 0.5], lw, '--');
add_local_coordinates(ax, Kp, -pi/6, len, 'black', lw, '-');
add_angle_indicator(ax, K, 0, pi/6, 1/4*len, '  \theta = \pi/6',...
    'black', [0, 0]);
add_angle_indicator(ax, Kp, 0, -pi/6, 1/4*len, '  \theta = -\pi/6',...
    'black', [0, 0]);

% add q vector symbol
thetaq = pi / 6 + pi / 4;
q_len = 0.5;
q = q_len * [cos(thetaq), sin(thetaq)];
quiver(ax, K(1), K(2), q(1), q(2), 'LineWidth', lw, 'Color', 'r');
add_angle_indicator(ax, K, 0, thetaq, 1/8*len, '  \theta_q = -5\pi/12',...
    'r', [-0.05, 0.2]);