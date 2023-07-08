tbg1 = TBG(cm=3, cn=2);
tbg2 = TBG(cm=7, cn=5);
tbg3 = TBG(cm=10, cn=7);
fig1 = tbg1.plot_real_space_lattice(0);
fig2 = tbg2.plot_real_space_lattice(0);
fig3 = tbg3.plot_real_space_lattice(0);

ax1 = fig1.Children;
ax2 = fig2.Children;
ax3 = fig3.Children;

axis(ax1, 'off');
axis(ax2, 'off');
axis(ax3, 'off');

tighten_margin(fig1, ax1);
tighten_margin(fig2, ax2);
tighten_margin(fig3, ax3);

% add circles to label orders
viscircles(ax1, [0, 0], norm(tbg1.a1l1), Color=[0.4660 0.6740 0.1880], LineStyle='-', LineWidth=3);
viscircles(ax1, 1/3*(2*tbg1.A2-tbg1.A1)', norm(tbg1.a1l1), Color=[0.4660 0.6740 0.1880], LineStyle='--', LineWidth=3);

viscircles(ax2, [0, 0], norm(tbg2.a1l1), Color=[0.4660 0.6740 0.1880], LineStyle='-', LineWidth=3);
viscircles(ax2, 1/2*tbg2.A2', norm(tbg2.a1l1), Color=[0.4940 0.1840 0.5560], LineStyle='-', LineWidth=3);

viscircles(ax3, [0, 0], norm(tbg3.a1l1), Color=[0.4660 0.6740 0.1880], LineStyle='-', LineWidth=3);
viscircles(ax3, 1/3*tbg3.A2', norm(tbg3.a1l1), Color=[0.4940 0.1840 0.5560], LineStyle='-', LineWidth=3);
viscircles(ax3, 2/3*tbg3.A2', norm(tbg3.a1l1), Color=[0.9290 0.6940 0.1250], LineStyle='-', LineWidth=3);

saveas(fig1, '..\..\figures\one_order.png');
saveas(fig2, '..\..\figures\two_order.png');
saveas(fig3, '..\..\figures\three_order.png');