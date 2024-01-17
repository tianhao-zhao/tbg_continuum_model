% single Dirac cone behavior
single_Bs = 1:0.01:20;
single_Es = zeros(9, length(single_Bs));
for i = 0:8
    cur_E = sqrt(2 * i * 1.05e-34 * (1e6)^2 * 1.6e-19 * single_Bs) / 1.6e-19;
    single_Es(i+1, :) = cur_E;
end

for run_ind = 1:15
    fig = figure();
    ax = axes(fig);
    load(sprintf('batch2run%d', run_ind));
    llrp = LLRP(eta=eta, vF=vF, N_ll=N_ll, theta=theta, w=w);
    Bs = sqrt(3) / 4 * llrp.ktheta^2 * llrp.hbar / (llrp.e * 2 * pi) * data(1, :);
    scatter(ax, Bs, data(2, :), 5, '*');
    hold(ax, 'on');    
    %plot(ax, single_Bs, single_Es, 'r');
    ylim(ax, [0, 0.6]);
    xlim(ax, [0,30]); 
    xlabel(ax, 'B(T)');
    ylabel(ax, 'E(eV)')
    ax_title = sprintf('theta=%d, w=%d, N_{ll}=%d, p_{max}=%d, q_{max}=%d', theta,...
        w/(1.6e-19 * 1e-3), N_ll, p_max, q_max);
    title(ax, ax_title);
end
