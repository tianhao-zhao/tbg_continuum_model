theta = 5;
vF = 1e6;
eta = 0.9;
cm = CM(theta=theta, vF=vF, eta=eta);
k_vanhove = -1/6 * (cm.gjs(:, 2) + cm.gjs(:, 3));
E = cm.calculate_e(k_vanhove);
disp('calculated van hove');
disp(E(length(E)/2)/cm.electron);
disp(E(length(E)/2 + 1) / cm.electron);
disp('predicted van hove');
alpha = cm.w / (cm.hbar * cm.vF * norm(cm.KD));
vF_renorm = (1 - 3*cm.eta*alpha^2)/(1 + (3 + 3*cm.eta)*alpha^2) * cm.vF;
disp(cm.hbar*vF_renorm*norm(k_vanhove)/(cm.electron));
disp('predicted with original vF')
disp(cm.hbar*cm.vF*norm(k_vanhove)/(cm.electron) - cm.w/cm.electron)

disp('predicted v renorm');
disp(vF_renorm);
disp('real v renorm');
E1 = cm.calculate_e([0, 0]);
E2 = cm.calculate_e(k_vanhove/100); 
deltaE = E2(length(E2)/2) - E1(length(E1) / 2);
true_renorm_vF = deltaE / (norm(k_vanhove/100) * cm.hbar);
disp(true_renorm_vF);