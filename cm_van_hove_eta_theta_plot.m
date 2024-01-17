vF = 1.1 * 1e6;
etas = 0:0.05:1;
colors = [];
for colorx = 0.3:0.7/length(thetas):1
    colors = [colors; [colorx, 0, 0]];
end

w = 80e-3 * 1.6e-19;
thetas = 1.5:0.1:2.5;
fig = figure();
ax = axes(fig);
fig2 = figure();
ax2 = axes(fig2);
for ii = 1:length(thetas)
    theta = thetas(ii);
    vhs = zeros(1, length(etas));
    vfres = zeros(1, length(etas));
    for ind = 1:length(etas)
        eta = etas(ind);
        cm = CM(theta=theta, vF=vF, eta=eta, w=w);
        k_vanhove = -1/6 * (cm.gjs(:, 2) + cm.gjs(:, 3));
        E = cm.calculate_e(k_vanhove);
        vhs(ind) = (abs(E(length(E)/2)) + abs(E(length(E)/2+1)))/2;
        vF_renorm = vhs(ind) / (cm.hbar * norm(k_vanhove));
        vfres(ind) = vF_renorm;
    end
    plot(ax, etas, vhs/(1.6e-19), Color=colors(ii, :));
    hold(ax, 'on');
    plot(ax2, etas, vfres, Color=colors(ii, :));
    hold(ax2, 'on');
end
