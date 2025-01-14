lw = 2;

tbg = TBG(cm=7,cn=6);
fig = tbg.plot_network;
ax = fig.Children(1);
hold(ax, 'on');

% plot g1 and g2
quiver(ax, [0, 0], [0, 0], [tbg.g1(1), tbg.g2(1)], [tbg.g1(2), tbg.g2(2)],...
    AutoScale='off', LineWidth=lw, Color='k', MaxHeadSize=1/sqrt(3));

% text g1 g2
text(ax, tbg.g1(1)/2, tbg.g1(2)/2, 'g_3 ', 'FontSize', 18,...
    'HorizontalAlignment','right');
text(ax, tbg.g2(1)/2, tbg.g2(2)/2, '   g_2', 'FontSize', 18);

% plot KD
quiver(ax, 0, 0, -tbg.KD(1), -tbg.KD(2), "off", LineWidth=lw, Color='k',...
    MaxHeadSize=1);
text(ax, -tbg.KD(1), -tbg.KD(2), sprintf(' K_D'), 'FontSize',18);

% plot t123
t_center = [0, 2*sqrt(3)*norm(tbg.KD)];
ts = norm(tbg.KD)*[[-1/2,-1/2,1];[sqrt(3)/2,-sqrt(3)/2, 0]];
quiver(ax, t_center(1)*[1,1,1], t_center(2)*[1,1,1], -ts(1,:), -ts(2,:),...
    "off",LineWidth=lw, Color='k',MaxHeadSize=1);
text(ax, -ts(1,1)+t_center(1), -ts(2,1)+t_center(2), ' -t_1', 'FontSize',18);
text(ax, -ts(1,2)+t_center(1), -ts(2,2)+t_center(2), ' -t_3', 'FontSize',18);
text(ax, -ts(1,3)+t_center(1), -ts(2,3)+t_center(2), '-t_2 ', 'FontSize',18,...
    'HorizontalAlignment','right');

axis(ax, 'off');



tighten_margin(fig, ax)
saveas(fig, '..\..\figures\g_network.png');