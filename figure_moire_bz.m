lw = 1.5;
a = 1;
fig = figure();
ax = axes(fig);

plot_single_hexagon(ax, [0, 0], a, 0, 0, lw, 'k');
KM = [-a, 0]';
KMp = [-a*cos(pi/3), -a*sin(pi/3)]';
M = (KM + KMp) / 2;
Gamma = [0, 0]';

l = [KM, Gamma, M, KMp];
line(ax, l(1, :), l(2, :), Color='k', LineStyle='--', LineWidth=lw);
text(ax, KM(1)+0.1*a, KM(2)+0.05*a, 'K_M', Color='k', FontSize=18*lw);
text(ax, KMp(1), KMp(2)+0.15*a, 'K''_M', Color='k', Fontsize=18*lw);
text(ax, M(1)-0.15*a, M(2)-0.15*a, 'M_M', Color='k', FontSize=18*lw);
text(ax, Gamma(1), Gamma(2), '\Gamma_M', FontSize=18*lw);

axis(ax, 'off');
tighten_margin(fig, ax);
saveas(fig, '..\..\figures\cm_mbz.png');