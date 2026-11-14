% Parameters
Tc = 32; % K
kB = 0.08617; % meV/K
c = physconst('LightSpeed') * 10^+2; % c in cm
T = [3 15 20 27]; % Passing T in K directly
kT = T * kB;
f_THz = linspace(0.09, 1, 300); % 1THz ~ 35cm^-1
hw = f_THz * 4.13567; % 1THz = 4.31567 meV
nT = numel(T);
nW = numel(hw);

% Parametri liberi modello
p = 0.675;                       % band weight 
gamma_1 = [139 148.1 151.3 145.4];                 % in cm^-1
gamma_2 = [147.8 149.2 150.1 147.9];                 % in cm^-1
Delta_1 = [3.282 2.512 2.317 1.561];                 % in meV
Delta_2 = [0.974 1.268 1.372 1.032];                 % in meV
sigma_tot = zeros(nT, nW); % due bande sigma_tot = p*sigma_1 + (1-p)sigma_2
% Iterate over T for BCS prediction
for i = 1:nT
    s_1 = zimmermann_finite_tau(hw, Delta_1(i), kT(i), gamma_1(i)*c*1e-12*4.13567); % all in meV
    s_2 = zimmermann_finite_tau(hw, Delta_2(i), kT(i), gamma_2(i)*c*1e-12*4.13567);
    sigma_tot(i, :) = p.*s_1 + (1 - p).*s_2;
end

% experimental data
sigma_dc = 0.65;               % non written , from the graph units?
%clean_wpd_csv('raw_sigma1_MgB2.csv', 'Real_cond_MgB2.csv', T)
%clean_wpd_csv('raw_sigma2_MgB2.csv', 'Img_cond_MgB2.csv', T)
s_real_data = readtable('Real_cond_MgB2.csv');
s_img_data = readtable('Img_cond_MgB2.csv');
%s_real_data.sigma = s_real_data.sigma/sigma_dc; % normalizzation
%s_img_data.sigma = s_img_data.sigma/sigma_dc;

% plot set up
x = 10^12*f_THz/c;
colors = lines(nT); 

% Plotting real part --------------------------------------
figure
hold on
for i = 1:nT
    col = colors(i, :); 
    % 1. BCS
    plot(x, real(sigma_tot(i, :)*sigma_dc), 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $T = %.1f$~K', T(i)));
    % 2. Data
    idx_exp = abs(s_real_data.T - T(i)) < 1e-6;
    if any(idx_exp)
        plot(s_real_data.nu(idx_exp), s_real_data.sigma(idx_exp), '--o', ...
            'Displayname',sprintf('Data $T = %.1f$~K', T(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Real part Conductivity of $MgB_2$ ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
ylim([-0.1 1.1])
xlim([3 35])
xlabel('Frequency (cm^{-1})'); 
%ylabel('\sigma_1/\sigma_n (\omega)');
ylabel('\sigma_2(10^{-3}(\Omega/\\sq)^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
%exportgraphics(gcf, 'BCS_MgB2_real_cond_finite_tau_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');


% Plotting imaginary part -----------------------------------
figure
hold on
for i = 1:nT
    col = colors(i, :);
    % 1. BCS
    plot(x, imag(sigma_tot(i, :)*sigma_dc), 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $T = %.1f$~K', T(i)));
    % 2. Data
    idx_exp = abs(s_img_data.T - T(i)) < 1e-6;
    if any(idx_exp)
        plot(s_img_data.nu(idx_exp), s_img_data.sigma(idx_exp), '--o', ...
            'DisplayName', sprintf('Data $T = %.1f$~K', T(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Imaginary part Conductivity of %MgB_2$ ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
ylim([-0.1 7])
xlim([3 35])
xlabel('Frequency (cm^{-1})'); 
%ylabel('\sigma_2/\sigma_n (\omega)');
ylabel('\sigma_2(10^{-3}(\Omega/sq)^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
%exportgraphics(gcf, 'BCS_MgB2_imm_cond_finite_tau_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');
