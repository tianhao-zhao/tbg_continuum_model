vFs = [0.866, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 1.3] * 1e6;
colors = [];
for colorx = 0.5:0.5/length(vFs):1
    colors = [colors; [colorx, 0, 0]];
end

w = 90e-3 * 1.6e-19;
thetas = 1.5:0.1:4;
eta = 0.6;
fig = figure();
ax = axes(fig);
fig2 = figure();
ax2 = axes(fig2);
for ii = 1:length(vFs)
    vF = vFs(ii);
    vhs = zeros(1, length(thetas));
    vfres = zeros(1, length(thetas));
    for ind = 1:length(thetas)
        theta = thetas(ind);
        cm = CM(theta=theta, vF=vF, eta=eta, w=w);
        k_vanhove = -1/6 * (cm.gjs(:, 2) + cm.gjs(:, 3));
        E = cm.calculate_e(k_vanhove);
        vhs(ind) = (abs(E(length(E)/2)) + abs(E(length(E)/2+1)))/2;
        vF_renorm = vhs(ind) / (cm.hbar * norm(k_vanhove));
        vfres(ind) = vF_renorm;
    end
    plot(ax, thetas, vhs/(1.6e-19), Color=colors(ii, :));
    hold(ax, 'on');
    plot(ax2, thetas, vfres, Color=colors(ii, :));
    hold(ax2, 'on');
end
