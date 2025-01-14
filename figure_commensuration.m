% from paper https://arxiv.org/ftp/arxiv/papers/2007/2007.03542.pdf
% m, n are integers

% % plot figures given a tbg object
% tbg = TBG(cm = 5, cn = 3);
% 
% % plot real space
% fig = figure();
% ax = axes(fig);
% plot_hexagons(ax, tbg.a1l1, tbg.a2l1, 'r', 2, 1, false, [0, 0], 1);
% plot_hexagons(ax, tbg.a1l2, tbg.a2l2, 'b', 2, 1, false, [0, 0], 1);
fig = figure();
ax = axes(fig);

a1 = [3/2, -sqrt(3)/2];
a2 = [3/2, sqrt(3)/2];

p_max = 30;
mn_lim = 30;
theta_aM = [];
for n = 1:mn_lim
    for m = (n+1):mn_lim
        theta_1 = acos((n^2 + 4*n*m + m^2)/(2*(n^2 + n*m + m^2)));
        aMl1 = n*a1 + m*a2;
        aM = norm(aMl1);
        theta_aM = [theta_aM;[theta_1*180/pi, aM/norm(a1)]];
        if (m==2 && n==1) || (m==3 && n==2) ||...
                (m==5 && n==3) || (m==4 && n==3)
            text(ax, theta_1*180/pi + 0.5, aM/norm(a1) -1, sprintf("%d, %d", m, n));
        end
        hold(ax, 'on');
    end
end

scatter(ax, theta_aM(:, 1), theta_aM(:, 2));
hold(ax, 'on');

theta = 0.01:0.001:1;
for p = 1:p_max
    plot(ax, theta*180/pi, p/2*1./sin(theta/2));
    hold(ax, 'on');
end

ax.XLim = [0, 30];
ax.YLim = [0, 55];
xlabel(ax, "Rotation Angle (Degree)", "FontSize", 18);
ylabel(ax, "|A_1|/|a_1|", "FontSize", 18);


box(ax, 'on');

tighten_margin(fig, ax);
saveas(fig, '..\..\figures\commensurate_angles.png');

