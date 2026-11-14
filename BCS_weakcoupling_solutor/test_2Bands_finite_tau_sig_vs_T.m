% Parameters
Tc = 32; % K
kB = 0.08617; % meV/K
c = physconst('LightSpeed') * 10^+2; % c in cm
T = linspace(3, 31.9, 200); % Passing T in K directly
kT = T * kB;
nu = [8.2 16];
hw = nu*1e-12*c * 4.13567; % 1THz = 4.31567 meV
nT = numel(T);
nW = numel(hw);

% Parametri liberi modello
p = 0.675;                       % band weight 
gamma_1 = 139;                 % in cm^-1
gamma_2 = 148;                 % in cm^-1
Delta0_1 = 3.282;                 % in meV
Delta0_2 = 0.974;                 % in meV
sigma_tot = zeros(nT, nW); % due bande sigma_tot = p*sigma_1 + (1-p)sigma_2
% Iterate over T for BCS prediction
for i = 1:nT
    Delta_1 = Delta0_1*tanh(1.74*sqrt(Tc/T(i) - 1));
    Delta_2 = Delta0_2*tanh(1.74*sqrt(Tc/T(i) - 1));
    s_1 = zimmermann_finite_tau(hw, Delta_1, kT(i), gamma_1*c*1e-12*4.13567); % all in meV
    s_2 = zimmermann_finite_tau(hw, Delta_2, kT(i), gamma_2*c*1e-12*4.13567);
    sigma_tot(i, :) = p.*s_1 + (1 - p).*s_2;
end

% experimental data
s_real_data = readtable('Real_MgB2_sig_vs_T.csv');
s_img_data = readtable('Img_MgB2_sig_vs_T.csv');
%s_real_data.sigma = s_real_data.sigma/sigma_dc; % normalizzation
%s_img_data.sigma = s_img_data.sigma/sigma_dc;

% plot set up
colors = lines(nW); 
error = 0.15;
% Plotting real part --------------------------------------
figure
hold on
for i = 1:nW
    col = colors(i, :); 
    % 1. BCS
    plot(T, real(sigma_tot(:, i)), 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $\\nu = %.1f$~$cm^{-1}$', nu(i)));
    % 2. Data
    idx_exp = abs(s_real_data.nu - nu(i)) < 1e-6;
    if any(idx_exp)
        plot(s_real_data.T(idx_exp), s_real_data.sigma(idx_exp), '--o', ...
            'Displayname',sprintf('Data $\\nu = %.1f$~~$cm^{-1}$', nu(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Real part Conductivity of MgB2 ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
%ylim([-0.1 1.5])
%xlim([3 35])
xlabel('T (K)'); 
ylabel('\sigma_1/\sigma_n (\omega)');
%ylabel('\sigma_2(10^{-3}(\Omega/sq)^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_2Bands_MgB2_real_cond_vs_T_fit.pdf', 'ContentType', 'vector');


% Plotting imaginary part -----------------------------------
figure
hold on
for i = 1:nW
    col = colors(i, :);
    % 1. BCS
    plot(T, imag(sigma_tot(:, i)), 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $\\nu = %.1f$~$cm^{-1}$', nu(i)));
    % 2. Data
    idx_exp = abs(s_real_data.nu - nu(i)) < 1e-6;
    if any(idx_exp)
        plot(s_img_data.T(idx_exp), s_img_data.sigma(idx_exp), '--o', ...
            'Displayname',sprintf('Data $\\nu = %.1f$~~$cm^{-1}$', nu(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Imaginary part Conductivity of MgB2 ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
%ylim([-0.1 7])
%xlim([3 35])
xlabel('T (K)');  
ylabel('\sigma_2/\sigma_n (\omega)');
%ylabel('\sigma_2(10^{-3}(\Omega/sq)^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_2Bands_MgB2_imm_cond_vs_T_fit.pdf', 'ContentType', 'vector');