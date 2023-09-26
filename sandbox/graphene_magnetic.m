% this code is to test 
% graphene landau level in magnetic field 
% doesn't depend on q or k

N = 30;
a = diag(sqrt(1:N-1), 1);
ad = diag(sqrt(1:N-1), -1);
EA = 0;
EB = 0;
hbar = 1.05e-34;
vF = 1e6;
e = 1.6e-19;
a_lattice = 1.4e-10;
K = 2*pi/a_lattice * 2*sqrt(3)/9;

B_max = 10;
% H = H(q)
E_q = zeros(2*N, B_max);
V_q = zeros(2*N, 2*N, B_max);
for b = 1:B_max
    lc = sqrt(hbar/(e*b));
    H_q = hbar*vF*[EA*eye(N), 1/lc*a; 1/lc*ad, EB*eye(N)];
    [V, E] = eig(H_q);
    E_q(:, b) = diag(E);
    V_q(:, :, b) = V;
end
fig_q = figure();
ax_q = axes(fig_q);

for i = 1:1:size(E_q, 1)
    plot(ax_q, 1:B_max, E_q(i,:));
    hold(ax_q, 'on');
end

% H = H(k)
% K_const = (10, 10)
Kx = 1/2*K;
Ky = sqrt(3)/2*K;
E_k = zeros(2*N, B_max);
V_k = zeros(2*N, 2*N, B_max);
for b = 1:B_max
    lc = sqrt(hbar/(e*b));
    H_k = hbar*vF*[EA*eye(N), 1/lc*a - (Kx - 1i*Ky)*eye(N);...
        1/lc*ad - (Kx + 1i*Ky)*eye(N), EB*eye(N)];
    [V, E] = eig(H_k);
    E_k(:, b) = diag(E);
    V_k(:, :, b) = V;
end
fig_k = figure();
ax_k = axes(fig_k);
for i = 1:size(E_k, 1)
    plot(ax_k, 1:B_max, E_k(i,:));
    hold(ax_k, 'on');
end


fig_com = figure();
ax_com = axes(fig_com);
for i = 1:20
    plot(ax_com, 1:B_max, E_q(i, :));
    hold(ax_com, 'on');
    plot(ax_com, 1:B_max, E_k(i, :));
    hold(ax_com, "on");
end

fig_com2= figure();
ax_com2 = axes(fig_com2);
for i = 1:5
    plot(ax_com2, 1:B_max, E_q(i+size(E_q, 1)/2, :));
    hold(ax_com2, 'on');
    plot(ax_com2, 1:B_max, E_k(i+size(E_q, 1)/2, :));
    hold(ax_com2, "on");
end

