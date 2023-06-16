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
theta_aM = zeros(20*20, 2);
for m = 1:40
    for n = m:40
        theta_1 = acos((n^2 + 4*n*m + m^2)/(2*(n^2 + n*m + m^2)));
        aMl1 = n*a1 + m*a2;
        aMl2 = m*a1 + n*a2;
        theta_2 = acos(aMl1*aMl2'/(norm(aMl1)*norm(aMl2)));
        aM = norm(aMl2);
        theta_aM(n+m*20, :) = [theta_1, aM];
        text(ax, theta_1, aM, sprintf("%d, %d", m, n));
        hold(ax, 'on');
    end
end

scatter(ax, theta_aM(:, 1), theta_aM(:, 2));
hold(ax, 'on');
p_max = 30;
theta = 0.05:0.001:1;
for p = 1:p_max
    plot(ax, theta, norm(a1)*p/2*1./sin(theta/2));
    hold(ax, 'on');
end
