% import data
fig1 = figure();
fig2 = figure();
ax2 = axes(fig2);
ax1 = axes(fig1);
scatter(ax1, data(1,:),data(2,:), 5, '*');
ylim(ax1, [-0.2, 0.2]);
xlim(ax1, [0, 1]);


% reconstruct a llrp to get magnetic field 
llrp = LLRP(eta=eta, vF=vF, N_ll=N_ll, theta=theta, w=w);
Bs = sqrt(3) / 4 * llrp.ktheta^2 * llrp.hbar / (llrp.e * 2 * pi) * data(1, :);

scatter(ax2, Bs, data(2, :), 5, '*');
hold(ax2, "on");
ylim(ax2, [-0, 0.3]);
xlim(ax2, [0, 20]); 

% single Dirac cone behavior
single_Bs = 1:0.01:20;
for i = 0:8
    cur_E = sqrt(2 * i * llrp.hbar * llrp.vF^2 * llrp.e * single_Bs) / llrp.e;
    plot(ax2, single_Bs, cur_E, 'r');
    hold(ax2, 'on');
end

ax_title = sprintf('theta=%d, w=%d, N_ll=%d, p_max=%d, q_max=%d', theta,...
    w/(1.6e-19 * 1e-3), N_ll, p_max, q_max);
title(ax2, ax_title);
title(ax1, ax_title);






