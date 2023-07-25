cm = 3;
cn = 2;
tbg = TBG(cm=cm, cn=cn);

% plot real space atoms
fig = tbg.plot_real_space_lattice();
ax = fig.CurrentAxes(1);

% unit cell
unit_cell = [[0;0], tbg.A1, tbg.A1+tbg.A2, tbg.A2, [0;0]];
line(ax, unit_cell(1, :), unit_cell(2, :), LineStyle='--', lineWidth=2, Color=[0.75, 0.75, 0.75]);
hold(ax, 'on');

% plot real space unit vector
quiver(ax, [0, 0], [0, 0], [tbg.a1l1(1), tbg.a2l1(1)], [tbg.a1l1(2), tbg.a2l1(2)],...
    AutoScale="off", Color='r', LineWidth=2);
quiver(ax, [0, 0], [0, 0], [tbg.a1l2(1), tbg.a2l2(1)], [tbg.a1l2(2), tbg.a2l2(2)],...
    AutoScale="off", Color='b', LineWidth=2);

quiver(ax, [0, 0], [0, 0], [tbg.A1(1), tbg.A2(1)], [tbg.A1(2), tbg.A2(2)],...
    'AutoScale','off', 'Color',[0.25, 0.25, 0.25], LineWidth=2);
hold(ax, 'on');

% text
text(ax, tbg.a1l1(1), tbg.a1l1(2)-0.25, ' a^{(1)}_1', Color='r', FontSize=24);
text(ax, tbg.a2l1(1), tbg.a2l1(2)-0.25, ' a^{(1)}_2', Color='r', FontSize=24);
text(ax, tbg.a1l2(1), tbg.a1l2(2)+0.25, ' a^{(2)}_1', Color='b', FontSize=24);
text(ax, tbg.a2l2(1), tbg.a2l2(2)+0.25, ' a^{(2)}_2', Color='b', FontSize=24);

text(ax, tbg.A1(1), tbg.A1(2), ' A_1', Color='k', FontSize=30);
text(ax, tbg.A2(1), tbg.A2(2), ' A_2', Color='k', FontSize=30);

% figure setup
ax.XLim = [-4.5, 13.5];
ax.YLim = [-2.5, 10.5];
axis(ax, 'off');
tighten_margin(fig, ax);

% save
saveas(fig, '..\..\figures\moire_real_space.png');