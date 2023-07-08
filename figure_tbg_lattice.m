% plot figures given a tbg object
cm = 5;
cn = 3;
tbg = TBG(cm=cm, cn=cn);

% plot real space
fig = figure();
ax = axes(fig);
plot_hexagons(ax, tbg.a1l1, tbg.a2l1, 'r', 2, 1, false, [0, 0], 1);
plot_hexagons(ax, tbg.a1l2, tbg.a2l2, 'b', 2, 1, false, [0, 0], 1);
% plot super unit cell
line(ax, [0, tbg.A1(1), tbg.A1(1)+tbg.A2(1), tbg.A2(1), 0],...
    [0, tbg.A1(2), tbg.A1(2) + tbg.A2(2), tbg.A2(2), 0], 'LineStyle', '--',...
    'Color', [0.5, 0.5, 0.5]);