fig = figure();
ax = axes(fig);

data_names = {'run2.mat', 'run4.mat', 'run5.mat'};
legends = {'40', '60', '20'};
for data_name = data_names
    load(string(data_name));
    llrp = LLRP(eta=eta, vF=vF, N_ll=N_ll, theta=theta, w=w);
    Bs = sqrt(3) / 4 * llrp.ktheta^2 * llrp.hbar / (llrp.e * 2 * pi) * data(1, :);
    scatter(ax, Bs, data(2, :), 5, '*');
    hold(ax, "on");
end
legend(legends);
ylim(ax, [0, 0.3]);
xlim(ax, [0,40]); 