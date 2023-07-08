close all;
addpath(genpath([pwd, '\distmesh']));
% note:
% all 1 dimen vectors are column vectors

% scirpt control params
% plot lattice space or not
display_mode = true;
% g network cutoff multiplier number
network_cutoff_number = 8;
% tolerance
tol = 1e-1;


% input variables
% k: defined in first MBZ in unrotated coordinates
% in this coordinate, gamma-M line is x axis
% layer1 is rotated down ctheta/2 (clockwise)
% layer2 is rotated up ctheta/2 (counterclockwise)
kx = 3;
ky = 2;
k = [kx; ky];

% physics constants
% a: distance between nearest A and B atoms
a = 1;
vf = 1;
hbar = 1;
% first order fourier component 
w = 0.0005;

% hopping matrix
w_phase = exp(1i*2*pi/3);
T1 = w * [[1, 1]; [1, 1]];
T2 = w * [[1, w_phase']; [w_phase, 1]];
T3 = w * [[1, w_phase]; [w_phase', 1]];
Ts = cat(3, T1, T2, T3);
% lattice parameters
% lattice unit vectors
a1 = a * [3/2; -sqrt(3)/2];
a2 = a * [3/2; sqrt(3)/2];
% displacement of nearest neighbours
a_delta1 = 1/3*(a1 + a2);
a_delta2 = a_delta1 - a1;
a_delta3 = a_delta1 - a2;
% reciprocal space unit vectors
b1 = 2 * pi / a * [1/3; -sqrt(3)/3];
b2 = 2 * pi / a * [1/3; sqrt(3)/3];
% K and K' (K prime) vectors
K = 2 * pi / a * [1/3; sqrt(3)/9];
Kp = 2 * pi / a * [1/3; -sqrt(3)/9];

% commensurate index and angle
cm = 3;
cn = 2;
ctheta = acos((cn^2 + 4*cn*cm + cm^2)/(2 * (cn^2 + cn*cm + cm^2)));
% rotation matrix
% note that layer2 is rotated counterclockwise ctheta
rm1 = [[cos(-ctheta/2), -sin(-ctheta/2)]; [sin(-ctheta/2), cos(-ctheta/2)]];
rm2 = [[cos(ctheta/2), -sin(ctheta/2)]; [sin(ctheta/2), cos(ctheta/2)]];
% rotated unit vectors
% a1l1: a1 vector in layer1
a1l1 = rm1 * a1;
a2l1 = rm1 * a2;
a1l2 = rm2 * a1;
a2l2 = rm2 * a2;
% a_delta1l1: displacement delta1 in layer1
a_delta1l1 = rm1 * a_delta1;
a_delta2l1 = rm1 * a_delta2;
a_delta3l1 = rm1 * a_delta3;

a_delta1l2 = rm2 * a_delta1;
a_delta2l2 = rm2 * a_delta2;
a_delta3l2 = rm2 * a_delta3;

% super lattice unit vectors
A1 = cn*a1l1 + cm*a2l1;
A2 = -cm*a1l1 + (cm+cn)*a2l1;
% super lattice reciprocal unit vectors
B1 = (cm+cn)/(cm^2+cn^2+cm*cn) * rm1*b1 + cm/(cm^2+cn^2+cm*cn) * rm1*b2;
B2 = -cm/(cm^2+cn^2+cm*cn) * rm1*b1 + cn/(cm^2+cn^2+cm*cn) * rm1*b2;
% plot real space lattice
if display_mode
    fig1 = figure(1);
    ax1 = axes(fig1);
    set(fig1, 'Name', sprintf('Real space (m,n) = (%d, %d)', cm, cn));
    plot_hexagons(ax1, a1l1, a2l1, 'r', 1, 1, false, highlight=1);
    plot_hexagons(ax1, a1l2, a2l2, 'b', 1, 1, false, highlight=1);
    plot_hexagons(ax1, A1, A2, [0.5, 0.5, 0.5], 1, 1, false);
end
% K and Kp in layer1 and 2
K1 = rm1 * K;
K1p = rm1 * Kp;
K2 = rm2 * K;
K2p = rm2 * Kp;
%display reciprocal lattice 
if display_mode
    fig2 = figure(2);
    ax2 = axes(fig2);
    set(fig2, 'Name', 'Reciprocal space');
    plot_hexagons(ax2, rm1*b1, -rm1*b2, 'r', 1, -1, true);
    plot_hexagons(ax2, rm2*b1, -rm2*b2, 'b', 1, -1, true);
    plot_hexagons(ax2, B1, -B2, [0.5, 0.5, 0.5], 3, 1, true);
end

% construct g vector network space
% unit vectors
KD = K1 - K2;
% Moire space unit vectors, they are just B1 and B2
g1 = [0; 2 * KD(2)];
g2 = [3 * KD(1); KD(2)];
% Moire space Ks and Ks' (Ksp)
Ks = [KD(1); -KD(2)];
Ksp = KD;

% search for all satisfied g_m vectors
% first find all valid g's
% displacement of q(k) to its three nearest neighbours
g_deltas = [KD, KD - g1, KD - g2];
% find all g vectors in cutoff range
search_limit = ceil(3 * network_cutoff_number);
[search_range_x, search_range_y] = meshgrid(-search_limit:search_limit);
g_network = [g1, g2] * [reshape(search_range_x, 1, []); reshape(search_range_y, 1, [])];
g_network = g_network(:, vecnorm(g_network + KD, 2, 1)<= network_cutoff_number * norm(KD));

if display_mode
    fig3 = figure(3);
    set(fig3, 'Name', 'Coupled states network');
    scatter(g_network(1, :), g_network(2, :), 'r', 'filled');
    hold on;
    scatter(g_network(1, :) + KD(1), g_network(2, :) + KD(2), 'b', 'filled');
    axis equal;
    quiver([0, 0, 0], [0, 0, 0], g_deltas(1, :), g_deltas(2, :), 'AutoScale', 'off');
    quiver([0, 0], [0, 0], [g1(1), g2(1)], [g1(2), g2(2)], 'AutoScale', 'off');
end

% TOSTUDY does the number of states in layer1 have to equal
% number of states in layer2?
 
M = size(g_network, 2);
% % hamiltonian matrix
% H = zeros(4*M);
% % construct intralayer hamiltonian
% % use kx and ky
% for ii = 1:M
%     cur_k = k + g_network(:, ii);
%     % in layer1
%     fk = exp(1i*cur_k'*a_delta1l1) + exp(1i*cur_k'*a_delta2l1) + ...
%         exp(1i*cur_k'*a_delta3l1);
%     % in layer2
%     cur_k_l2 = cur_k + KD;
%     fk_l2 = exp(1i*cur_k_l2'*a_delta1l2) + exp(1i*cur_k_l2'*a_delta2l2) + ...
%         exp(1i*cur_k_l2'*a_delta3l2);
%     H(2*ii - 1, 2*ii) = fk;
%     H(2*ii, 2*ii - 1) = fk';
%     H(2*ii - 1 + 2*M, 2*ii + 2*M) = fk_l2;
%     H(2*ii + 2*M, 2*ii - 1 + 2*M) = fk_l2';
% end
% % construct interlayer hamiltonian
% for ii = 1:M
%     for jj = 1:M
%         % ii in layer1, jj in layer2
%         dist = (g_network(:, jj) + KD) - g_network(:, ii);
%         for dd = 1:3
%             if norm(g_deltas(:, dd) - dist) <= tol
%                 % not sure if conjugate should be reversed 
%                 H(ii*2-1:ii*2, jj*2-1+2*M:jj*2+2*M) = Ts(:, :, dd)';
%                 H(jj*2-1+2*M:jj*2+2*M, ii*2-1:ii*2) = Ts(:, :, dd);
%                 % for test if all connections are found
%                 figure(fig3);
%                 cur_line = [g_network(:, ii), g_network(:, jj) + KD];
%                 line(cur_line(1, :), cur_line(2, :));
%             end
%         end
%     end
% end
% 
% % solve for eigenvalues for one momentum (q)
% E = eig(H);

% solve for all q's and connect
 
% plot dispersion in hexagonsal region
figure(5)
[k_space, t] = distmesh2d(@(p) sqrt(sum(p.^2, 2))-1.5*norm(Ks), @huniform,...
    norm(Ks)/30, [-2*norm(Ks), -2*norm(Ks); 2*norm(Ks), 2*norm(Ks)], []);

% plot dispersion for K gamma M K lines
K_origin = [0; 0];
M_point = (Ks + Ksp) / 2;
% when origin is K point, gamma is Ksp
gamma = Ksp;
grain = (0:100) ./ 100;
k_line = [(gamma - K_origin) * grain + K_origin, ...
    (M_point - gamma) * grain + gamma, ...
    (K_origin - M_point) * grain + M_point];
k_line_length = [vecnorm((gamma - K_origin) * grain, 2), ...
    vecnorm((M_point - gamma) * grain, 2) + norm(gamma - K_origin), ...
    vecnorm((K_origin - M_point) * grain, 2) + norm(gamma - K_origin) + norm(M_point - gamma)];



%test sandbox
% plot K gamma M K line
% origin is K point
tic
k_line = k_line';
H_inter = zeros(2*M);
for ii = 1:M
    for jj = 1:M
        % ii in layer1, jj in layer2
        dist = (g_network(:, jj) + KD) - g_network(:, ii);
        for dd = 1:3
            if norm(g_deltas(:, dd) - dist) <= tol
                % not sure if conjugate should be reversed 
                H_inter(ii*2-1:ii*2, jj*2-1:jj*2) = Ts(:, :, dd)';
                % for test if all connections are found
                % figure(fig3);
                % cur_line = [g_network(:, ii), g_network(:, jj) + KD];
                % line(cur_line(1, :), cur_line(2, :));
            end
        end
    end
end

all_E = zeros(4*M, size(k_line, 1));
for index = 1:size(k_line, 1)
    k = k_line(index, :)';
    M = size(g_network, 2);
    % hamiltonian matrix
    H_intra_1 = zeros(2*M);
    H_intra_2 = zeros(2*M);
    % construct intralayer hamiltonian
    % use kx and ky
    for ii = 1:M
        cur_k_l1 = k + g_network(:, ii);
        cur_k_l2 = cur_k_l1 + KD;
        % % in layer1
        % fk = exp(1i*cur_k'*a_delta1l1) + exp(1i*cur_k'*a_delta2l1) + ...
        %     exp(1i*cur_k'*a_delta3l1);
        % % in layer2
        % % we are using k (origin at gamma point, no need to add KD)
        % cur_k_l2 = cur_k;
        % fk_l2 = exp(1i*cur_k_l2'*a_delta1l2) + exp(1i*cur_k_l2'*a_delta2l2) + ...
        %     exp(1i*cur_k_l2'*a_delta3l2);

        H_intra_1(2*ii - 1, 2*ii) = hbar * vf * ...
            (cur_k_l1(1) - 1i * cur_k_l1(2)) * exp(-1i * ctheta/2);
        H_intra_1(2*ii, 2*ii - 1) = hbar * vf * ...
            (cur_k_l1(1) + 1i * cur_k_l1(2)) * exp(1i * ctheta/2);
        H_intra_2(2*ii - 1, 2*ii) = hbar * vf * ...
            (cur_k_l2(1) - 1i * cur_k_l2(2)) * exp(1i * ctheta/2);
        H_intra_2(2*ii, 2*ii - 1) = hbar * vf * ...
            (cur_k_l2(1) + 1i * cur_k_l2(2)) * exp(-1i * ctheta/2);
    end
    H = [H_intra_1, H_inter; H_inter', H_intra_2];
    % solve for eigenvalues for one momentum (q)
    E = eig(H);
    all_E(:,index) = E;
end
% plot 
figure(7)
for ii = 1:4*M
    plot(k_line_length, all_E(ii, :))
    hold on
end
vf_e = (all_E(2*M+1, 5) - all_E(2*M+1, 1)) / ...
    (k_line_length(5) - k_line_length(1));
fprintf('effective velocity = %d', vf_e) 
toc



% plot 3d dispersion
% tic
% % interlayer H block, up right block
% H_inter = zeros(2*M);
% for ii = 1:M
%     for jj = 1:M
%         % ii in layer1, jj in layer2
%         dist = (g_network(:, jj) + KD) - g_network(:, ii);
%         for dd = 1:3
%             if norm(g_deltas(:, dd) - dist) <= tol
%                 % not sure if conjugate should be reversed 
%                 H_inter(ii*2-1:ii*2, jj*2-1:jj*2) = Ts(:, :, dd)';
%                 % for test if all connections are found
%                 % figure(fig3);
%                 % cur_line = [g_network(:, ii), g_network(:, jj) + KD];
%                 % line(cur_line(1, :), cur_line(2, :));
%             end
%         end
%     end
% end
% 
% % from now on, the k is q, q = k - K
% all_E = zeros(4*M, size(k_space, 1));
% for index = 1:size(k_space, 1)
%     k = k_space(index, :)';
%     M = size(g_network, 2);
%     % hamiltonian matrix
%     H_intra_1 = zeros(2*M);
%     H_intra_2 = zeros(2*M);
%     % construct intralayer hamiltonian
%     % use kx and ky
%     for ii = 1:M
%         cur_k_l1 = k + g_network(:, ii);
%         cur_k_l2 = cur_k_l1 + KD;
%         % % in layer1
%         % fk = exp(1i*cur_k'*a_delta1l1) + exp(1i*cur_k'*a_delta2l1) + ...
%         %     exp(1i*cur_k'*a_delta3l1);
%         % % in layer2
%         % % we are using k (origin at gamma point, no need to add KD)
%         % cur_k_l2 = cur_k;
%         % fk_l2 = exp(1i*cur_k_l2'*a_delta1l2) + exp(1i*cur_k_l2'*a_delta2l2) + ...
%         %     exp(1i*cur_k_l2'*a_delta3l2);
% 
%         H_intra_1(2*ii - 1, 2*ii) = hbar * vf * ...
%             (cur_k_l1(1) - 1i * cur_k_l1(2)) * exp(-1i * ctheta/2);
%         H_intra_1(2*ii, 2*ii - 1) = hbar * vf * ...
%             (cur_k_l1(1) + 1i * cur_k_l1(2)) * exp(1i * ctheta/2);
%         H_intra_2(2*ii - 1, 2*ii) = hbar * vf * ...
%             (cur_k_l2(1) - 1i * cur_k_l2(2)) * exp(1i * ctheta/2);
%         H_intra_2(2*ii, 2*ii - 1) = hbar * vf * ...
%             (cur_k_l2(1) + 1i * cur_k_l2(2)) * exp(-1i * ctheta/2);
%     end
%     H = [H_intra_1, H_inter; H_inter', H_intra_2];
%     % solve for eigenvalues for one momentum (q)
%     E = eig(H);
%     all_E(:,index) = E;
% end
% toc
% % figure(4)
% % for ii = 1:size(all_E, 1)
% %     surf(k_space(:, 1), k_space(:, 2), all_E(ii, :))
% %     hold on
% % end
% 
% figure(8)
% trisurf(t, k_space(:, 1), k_space(:, 2), all_E(2*M, :))
% hold on
% trisurf(t, k_space(:, 1), k_space(:, 2), all_E(2*M+1, :))


% m = 40;
% M = 1;
% all_E = zeros(2*M, (2*m+1)*(2*m+1));
% for kxx = -m:m
%     for kyy = -m:m
%         kx = kxx * norm(B1) / 20;
%         ky = kyy * norm(B1) / 20;
%         H = zeros(2*M);
%         k = [kx; ky];
%         % construct intralayer hamiltonian
%         % use kx and ky
%         for ii = 1:M
%             cur_k = k;
%             % % in layer1
%             % fk = exp(i*cur_k'*a_delta1l1) + exp(i*cur_k'*a_delta2l1) + ...
%             %     exp(i*cur_k'*a_delta3l1);
%             % in layer2
%             cur_k_l2 = cur_k;
%             fk_l2 = exp(1i*cur_k_l2'*a_delta1l2) + exp(1i*cur_k_l2'*a_delta2l2) + ...
%                 exp(1i*cur_k_l2'*a_delta3l2);
%             % H(2*ii - 1, 2*ii) = fk;
%             % H(2*ii, 2*ii - 1) = fk';
%             H(2*ii - 1, 2*ii) = fk_l2;
%             H(2*ii, 2*ii - 1) = fk_l2';
%         end
%         E = eig(H);
%         all_E(:,(kxx+m)*(2*m+1)+(kyy+m)+1) = E;
%     end
% end
% figure(6)
% for ii = 1:size(all_E, 1)
%     surf(-m:m, -m:m, reshape(all_E(ii, :), 2*m+1, 2*m+1))
%     hold on
% end
