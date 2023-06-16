% plot figures given a tbg object
tbg = TBG(cm = 3, cn = 2);

% plot real space
fig = figure();
ax = axes(fig);
plot_hexagons(ax, tbg.a1l1, tbg.a2l1, 'r', 2, 1, false, [0, 0], 1);
plot_hexagons(ax, tbg.a1l2, tbg.a2l2, 'b', 2, 1, false, [0, 0], 1);