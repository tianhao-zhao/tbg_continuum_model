theta = 5;
vF = 1e6;
cc = 1;
cm = CM(theta=theta, vF=vF, cc = 1);
k_vanhove = -1/6 * (cm.gjs(:, 2) + cm.gjs(:, 3));
E = cm.calculate_e(k_vanhove);
disp('calculated van hove');
disp(E(length(E)/2)/cm.electron);
disp('predicted van hove');
alpha = cm.w / (cm.hbar * cm.vF * norm(cm.KD));
vF_renorm = (1 - 3*cm.cc*alpha^2)/(1 + (3 + 3*cm.cc)*alpha^2) * cm.vF;
disp(cm.hbar*vF_renorm*norm(cm.KD)/(2*cm.electron) - cm.w/cm.electron);
disp('predicted with original vF')
disp(cm.hbar*cm.vF*norm(cm.KD)/(2*cm.electron) - cm.w/cm.electron)
