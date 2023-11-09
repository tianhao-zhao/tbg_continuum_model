thetas = [10, 5, 2.5, 1.1];
cm = CM();

fig = figure();
tlo = tiledlayout(fig, 3, length(thetas));

% vF = 0.8e6, cc = 1
for i = 1:length(thetas)
    theta = thetas(i);
    cm.init(theta=theta);
    [~, ax] = cm.plot_kgammamk();
    ax.Parent = tlo;
    ax.Layout.Tile = i;
    cur_title = sprintf('\\theta = %0.1f deg, v_F = %0.2e m/s, \\eta = %0.1f',...
        cm.theta, cm.vF, cm.cc);
    subtitle(ax, cur_title);
end

% vF = 0.8e6, cc = 0.6
cc = 0.6;
for i = 1:length(thetas)
    theta = thetas(i);
    cm.init(theta=theta, cc=cc);
    [~, ax] = cm.plot_kgammamk();
    ax.Parent = tlo;
    ax.Layout.Tile = length(thetas)+i;
    cur_title = sprintf('\\theta = %0.1f deg, v_F = %0.2e m/s, \\eta = %0.1f',...
        cm.theta, cm.vF, cm.cc);
    subtitle(ax, cur_title);
end

% vF = 1.05e6, cc = 0.6
cc = 0.6;
vF = 1.05e6;
for i = 1:length(thetas)
    theta = thetas(i);
    cm.init(theta=theta, cc=cc, vF=vF);
    [~, ax] = cm.plot_kgammamk();
    ax.Parent = tlo;
    ax.Layout.Tile = 2*length(thetas)+i;
    cur_title = sprintf('\\theta = %0.1f deg, v_F = %0.2e m/s, \\eta = %0.1f',...
        cm.theta, cm.vF, cm.cc);
    subtitle(ax, cur_title);
end

ylabel(tlo, 'Energy (eV)');
tlo.TileSpacing = "tight";
tlo.Padding = 'tight';

set(fig, 'Position', [1, 1, 1000, 2000]);
saveas(fig, '..\..\figures\cm_kgammamkp.png');