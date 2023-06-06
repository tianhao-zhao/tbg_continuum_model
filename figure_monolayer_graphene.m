% monolayer graphene fig
% 1. real space
% 2. reciprocal space
% 3. 3d band structure
a = 1;
% base linewidth
lw = 1.5;

fig1 = figure();
fig2 = figure();
fig3 = figure();
% in pixels
fig1.Position = [0, 0, 1500, 1000];
fig2.Position = [0, 0, 1000, 1000];
fig3.Position = [0, 0, 1500, 1000];

ax1 = axes(fig1);
ax2 = axes(fig2);
ax3 = axes(fig3);

%ax1.Position = [0, 0, 1, 1];

% real space
a1 = a * [3/2; sqrt(3)/2];
a2 = a * [3/2; -sqrt(3)/2];
as = [a1, a2];
center1 = a * [-1; 0];
[xx, yy] = meshgrid(-3:3);
xx = xx(:);
yy = yy(:);
xy = [xx, yy]';
centers = [a1, a2] * xy + center1;

for ii = 1:size(centers, 2)
    plot_single_hexagon(ax1, centers(:, ii), a, 0, 1, lw);
end

% add a1 and a2 vectors
quiver(ax1, [0, 0], [0, 0], as(1, :), as(2, :), AutoScale='off', LineWidth=2*lw);
hold(ax1, 'on');
quiver(ax1, [0, 0, 0], -sqrt(3)*[1, 1, 1],...
    [1, -1/2, -1/2], [0, sqrt(3)/2, -sqrt(3)/2], AutoScale='off', LineWidth=2*lw);
hold(ax1, 'on');

text(ax1, -0.4, 0, 'A', FontSize=18*lw);
hold(ax1, 'on');
text(ax1, a + 0.2, 0, 'B', FontSize=18*lw);
hold(ax1, 'on');
text(ax1, 0.5, 0.7, 'a_2', FontSize=18*lw);
hold(ax1, 'on');
text(ax1, 0.5, -0.5, 'a_1', FontSize=18*lw);
text(ax1, 0.5, -1.5, '\delta_1', FontSize=18*lw);
text(ax1, -0.25, -1.2, '\delta_2', FontSize=18*lw);
text(ax1, -0.25, -2.3, '\delta_3', FontSize=18*lw);
hold(ax1, 'off');

axis(ax1, 'off');
ax1.XLim = [-3.1, 3.2];
ax1.YLim = [-4.1, 2.2];

%ax1.Title.String = '(a) Crystal lattice';

% reciprocal space
plot_single_hexagon(ax2, [0, 0], a, pi/6, 0, lw);
b1 = a * sqrt(3)/2 * [1; sqrt(3)];
b2 = a * sqrt(3)/2 * [1; -sqrt(3)];
bs = [b1, b2];

quiver(ax2, [0, 0], [0, 0], bs(1, :), bs(2, :), 'AutoScale', 'off', LineWidth=2*lw);
hold(ax2, 'on');

text(ax2, 0.4, 1.1, 'b_1', FontSize=18*lw);
text(ax2, 0.4, -1, 'b_2', FontSize=18*lw);
text(ax2, 0.9, 0.5, 'K', FontSize=18*lw);
text(ax2, 0.9, -0.5, "K'", FontSize=18*lw);
hold(ax2, 'off');

axis(ax2, 'off');
%ax2.PlotBoxAspectRatio = [0.5, 1, 1];
%ax2.Title.String = '(b) First Brillouin zone in reciprocal space';
% 3d band
kmax = 2.5;
gran = 100;
nmax = kmax * gran;
d1 = [1; 0];
d2 = [-1/2; sqrt(3)/2];
d3 = [-1/2; -sqrt(3)/2];
[kx, ky] = meshgrid(-nmax:nmax);
kx = kx(:) ./ gran;
ky = ky(:) ./ gran;
Ep = zeros(size(kx, 1), 1);
for ii = 1:size(kx, 1)
    ck = [kx(ii), ky(ii)];
    Ep(ii, 1) = abs(exp(1i*ck*d1) + exp(1i*ck*d2) +...
        exp(1i*ck*d3));
end
Ep = reshape(Ep, 2*nmax + 1, 2*nmax + 1);
xx = (-nmax:nmax) ./gran;
yy = xx;
surf(ax3, xx, yy, Ep, FaceColor='interp', LineStyle='none');
hold(ax3, 'on');
surf(ax3, xx, yy, -Ep, FaceColor="interp", LineStyle="none");
hold(ax3, 'off');
axis(ax3, 'equal');
%ax3.Title.String = '(c) Band structure';
ax3.XLabel.String = 'k_x';
ax3.XLabel.FontSize = 18*lw;
ax3.YLabel.String = 'k_y';
ax3.YLabel.FontSize = 18*lw;

% save to figures/thesis_figures
addpath(genpath([pwd, '\figures']));
saveas(fig1, '..\..\figures\monolayer_crystal_structure.png');
saveas(fig2, '..\..\figures\monolayer_reciprocal_lattice.png');
saveas(fig3, '..\..\figures\monolayer_bands.png');
