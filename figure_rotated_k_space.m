lw = 1.5;
ctheta = 8 / 180 * pi;
a = 1;
rm1 = [[cos(-ctheta), -sin(-ctheta)]; [sin(-ctheta), cos(-ctheta)]];
rm2 = [[cos(ctheta), -sin(ctheta)]; [sin(ctheta), cos(ctheta)]];

fig = figure();
ax = axes(fig);

% reciprocal space
plot_single_hexagon(ax, [0, 0], a, pi/6 + ctheta, 0, lw, 'b');
b1 = rm1* (a * sqrt(3)/2 * [1; -sqrt(3)]);
b2 = rm1* (a * sqrt(3)/2 * [1; sqrt(3)]);
bs = [b1, b2];

plot_single_hexagon(ax, [0, 0], a, pi/6 -ctheta, 0, lw, 'r');
b1_2 = rm2* (a * sqrt(3)/2 * [1; -sqrt(3)]);
b2_2 = rm2* (a * sqrt(3)/2 * [1; sqrt(3)]);
bs_2 = [b1_2, b2_2];

% b^1 and b^2 vectors
quiver(ax, [0, 0], [0, 0], bs(1, :), bs(2, :), 'r', 'AutoScale', 'off', LineWidth=2*lw);
hold(ax, 'on');
quiver(ax, [0, 0], [0, 0], bs_2(1, :), bs_2(2, :), 'b', 'AutoScale', 'off', LineWidth=2*lw);
hold(ax, 'on');

% calculate and plot mbz 
K1 = rm1 * (a * [sqrt(3)/2; 1/2]);
K2 = rm2 * (a * [sqrt(3)/2; 1/2]);
a_mbz = norm(K2 - K1);
mbz_center = K1 + a_mbz * [1/2; sqrt(3)/2];
plot_single_hexagon(ax, mbz_center, a_mbz, 0, 0, lw, 'k');

% Kd vector
% KD = K2 - K1;
% quiver(ax, K1(1), K1(2), KD(1), KD(2), 'k', 'AutoScale', 'off', LineWidth=2*lw);

% text for b vectors
text(ax, b1(1), b1(2)+0.1, ' b_1^{(1)}', Color='r', FontSize=18*lw);
text(ax, b2(1), b2(2), ' b_2^{(1)}', Color='r', FontSize=18*lw);
text(ax, b1_2(1), b1_2(2), ' b_1^{(2)}', Color='b', FontSize=18*lw);
text(ax, b2_2(1), b2_2(2)-0.1, ' b_2^{(2)}', Color='b', FontSize=18*lw);

% text for K points
text(ax, K1(1), K1(2), ' K^{(1)}', Color='r', FontSize=18*lw);
text(ax, K2(1), K2(2), ' K^{(2)}', Color='b', FontSize=18*lw);

% text for three b's
text(ax, K1(1)-b2(1), K1(2)-b2(2), sprintf('K^{(1)}-b_2^{(1)}\n=K^{(1)}-G^{(1)}_2'),...
    Color=[0.5, 0.5, 0.5],...
    FontSize=18*lw, HorizontalAlignment='right');
text(ax, K1(1)-b1(1)-b2(1), K1(2)-b1(2)-b2(2), sprintf('K^{(1)}-b^{(1)}_1-b^{(1)}_2\n=K^{(1)}-G^{(1)}_3'),...
    Color=[0.5, 0.5, 0.5], FontSize=18*lw, HorizontalAlignment='right');
text(ax, K1(1), K1(2)-0.1, sprintf('\n K^{(1)}-0=K^{(1)}-G^{(1)}_1'), Color=[0.5, 0.5, 0.5], FontSize=18*lw);

% arrow for three b's
quiver(ax, [K1(1), K1(1)], [K1(2), K1(2)], [-b2(1), -b2(1)-b1(1)], [-b2(2), -b2(2)-b1(2)],...
    AutoScale='off', Color=[0.5, 0.5, 0.5], LineWidth=2*lw);

% arrow for g2 and g3
g2 = b2 - b2_2;
g3 = (b1+b2) - (b1_2+b2_2);
quiver(ax, [mbz_center(1), mbz_center(1)], [mbz_center(2), mbz_center(2)],...
    [g2(1), g3(1)], [g2(2), g3(2)], AutoScale='off', Color='k', LineWidth=2*lw);
text(ax, mbz_center(1)+g2(1)/2, mbz_center(2)+g2(2)/2, '   g_2', 'FontSize',18*lw);
text(ax, mbz_center(1)+g3(1), mbz_center(2)+g3(2), sprintf('\ng_3'), 'FontSize',18*lw);

hold(ax, 'off');
axis(ax, 'off');
tighten_margin(fig, ax);

saveas(fig, '..\..\figures\cm_rotated_k_space.png');



