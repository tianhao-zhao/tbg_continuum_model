tic;
w = 1 * 110e-3 * 1.6e-19;
theta = 2.5;
vF = 1e6;
eta = 0.9;
N_ll = 15;
p_max = 10;
q_max = 10;
p_range = 1:1:p_max;
q_range = 1:1:q_max;
coprime_pq_pairs = [];

for p = p_range
    for q = q_range
        if iscoprime([p, q]) && p >= q
            coprime_pq_pairs = [coprime_pq_pairs, [p; q]];
        end
    end
end 

fig1 = figure();              
ax1 = axes(fig1);
data = [];

parfor ind = 1:length(coprime_pq_pairs)
    p = coprime_pq_pairs(1, ind);
    q = coprime_pq_pairs(2, ind);
    llrp = LLRP(theta=theta, N_ll=N_ll, vF=vF, eta=eta, w=w);
    llrp.precompute_init(p, q);
    llrp.construct_H(0, 0);
    cur_E = eig(llrp.H)' / llrp.e;
    ratio = q / p;
    cur_B = llrp.B;
    data = [data, [repelem(ratio, length(cur_E)); cur_E]];
end
scatter(data(1,:), data(2,:), 5, 'k', 'filled');
toc;