% fig1 = figure();
% ax1 = axes(fig1);
% plot(ax1, [1:100], [1:100]*2);
% 
% fig2 = figure();
% ax2 = axes(fig2);
% plot(ax2, [1:100], [1:100]*10);
% 
% ax_array = [ax1, ax2];
% 
% fig = figure();
% tlo = tiledlayout(2, 2);
% for i = 1:2
%     ax_array(i).Parent = tlo;
%     ax_array(i).Layout.Tile = i*2;
%     cur_title = sprintf('\\theta = %0.1f {\circ}, v_F = %0.2e m/s, \\alpha = %0.1f',...
%         1,2,3);
%     subtitle(ax_array(i), cur_title);
% end

% H = zeros(4, 4);
% ticBytes(gcp);
% parfor i = 1:4
%     for j = 1:4
%         H(i, j) = 1;
%     end
% end
% tocBytes(gcp);

tic;
foo = 1.6;
dp = zeros(20, 20);
for n =1:20
    for np = 1:n
        dp(n, np) = laguerreL(np, n - np, foo);
    end
end
toc;

