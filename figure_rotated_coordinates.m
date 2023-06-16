% linewidth
lw = 1.5;
len = 1;
font_scale = 18;

% K K' points
K = [sqrt(3)/2, 1/2];
Kp = [sqrt(3)/2, -1/2];

fig = figure();
ax = axes(fig);
plot_single_hexagon(ax, [0, 0], 1, pi/6, 0, lw);
add_local_coordinates(ax, K, 0, len, [0.5, 0.5, 0.5], lw, '--');
add_local_coordinates(ax, K, pi/6, len, 'red', lw, '-');
add_local_coordinates(ax, Kp, 0, len, [0.5, 0.5, 0.5], lw, '--');
add_local_coordinates(ax, Kp, -pi/6, len, 'blue', lw, '-');
add_angle_indicator(ax, K, 0, pi/6, 1/4*len, '  \theta = \pi/6',...
    [0.5, 0.5, 0.5], [0, 0], lw/2);
add_angle_indicator(ax, Kp, 0, -pi/6, 1/4*len, '  \theta = -\pi/6',...
    [0.5, 0.5, 0.5], [0, 0], lw/2);

% add q vector symbol
thetaq = pi / 6 + pi / 4;
q_len = 0.5;
q = q_len * [cos(thetaq), sin(thetaq)];
quiver(ax, K(1), K(2), q(1), q(2), 'LineWidth', lw, 'Color', 'r',...
    'AutoScale', 'off');
add_angle_indicator(ax, K, pi/6, thetaq, 1/8*len, '  \theta_q - \theta',...
    'r', [-0.05, 0.05], lw/2);
text(ax, q(1) + K(1), q(2) + K(2), ' P', 'Color', 'black',...
    'FontSize', font_scale*lw);


% add q projection
line(ax, [q(1) + K(1), q_len * cos(pi/4) * cos(pi/6) + K(1)],...
    [q(2) + K(2), q_len * cos(pi/4) * sin(pi/6) + K(2)],...
    'LineWidth', lw, 'Color', 'r', 'LineStyle', '--');
line(ax, [q(1) + K(1), -q_len * sin(pi/4) * sin(pi/6) + K(1)],...
    [q(2) + K(2), q_len * sin(pi/4) * cos(pi/6) + K(2)],...
    'LineWidth', lw, 'Color', 'r', 'LineStyle', '--')

% add center axis
add_local_coordinates(ax, [0, 0], 0, 3, 'black', lw, '-');
text(ax, K(1) - 0.8, K(2), 'Global', 'Color', [0.5, 0.5, 0.5],...
    'FontSize', font_scale*lw);
text(ax, K(1) - 0.8, K(2) - 0.3, 'Local', 'Color', 'r',...
    'FontSize', font_scale*lw);

tighten_margin(fig, ax);

axis(ax, 'off');
saveas(fig, '..\..\figures\rotating_coordinates.png');