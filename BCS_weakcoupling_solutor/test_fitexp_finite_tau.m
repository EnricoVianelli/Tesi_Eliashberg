% Parameters
Tc     = 8.3;                           % K
kB     = 0.08617;                       % meV/K
alpha = 4.1;                           % coupling 2*D0/KbTc
Delta0 = 0.5*alpha*kB*Tc;               % meV, BCS weak coupling
c = physconst('LightSpeed') * 10^+2;    % c in cm
gamma = 150;                            % in cm^-1
T = [4.5 7];                          % Passing T in K directly
kT = T * kB;

f_THz = linspace(0.01, 2, 400);   % 1THz ~ 35cm^-1
hw = f_THz * 4.13567;
x = 10^12*f_THz/c;
nT = numel(T);
nW = numel(hw);

sigma1 = zeros(nT, nW);
sigma2 = zeros(nT, nW);
sigma1_gamma = zeros(nT, nW);
sigma2_gamma = zeros(nT, nW);

% Iterate over T for BCS prediction
for i = 1:nT
    Delta = Delta0*tanh(1.74*sqrt(Tc/T(i) - 1));
    % no gamma
    [s1, s2] = mattis_bardeen(hw, Delta, kT(i));
    sigma1(i, :) = s1; sigma2(i, :) = s2;
    % gamma finito
    [s1_gamma, s2_gamma] = mb_finite_tau(hw, Delta, kT(i), gamma*c*1e-12*4.13567);
    sigma1_gamma(i, :) = s1_gamma; sigma2_gamma(i, :) = s2_gamma;
end

% experimental data
sigma_dc = 2.4; % units sigma [10^5 ohm^-1 cm^-1] at T = 9 K
raw = readtable('nb_digitized_raw.csv');
s1_data = raw(raw.panel==1,:);  s2_data = raw(raw.panel==2,:);
%s1_data.sigma = s1_data.sigma/sigma_dc;           % normalizzation
%s2_data.sigma = s2_data.sigma/sigma_dc;

% colors
colors = lines(nT); 

% Plotting real part --------------------------------------
figure
hold on
for i = 1:nT
    col = colors(i, :); 
    % 1. BCS
    plot(x, sigma1(i, :)*sigma_dc, 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $T = %.1f$~K', T(i)));
    % 2. BCS with gamma
    plot(x, sigma1_gamma(i, :)*sigma_dc, '--', 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS-gamma $T = %.1f$~K', T(i)));
    % 3. Data
    idx_exp = abs(s1_data.T - T(i)) < 1e-6;
    if any(idx_exp)
        plot(s1_data.nu(idx_exp), s1_data.sigma(idx_exp), 'o', ...
            'Displayname',sprintf('Data $T = %.1f$~K', T(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Real part Conductivity of Nb ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
ylim([-0.1 3])
xlabel('Frequency (cm^{-1})'); 
%ylabel('\sigma_1/\sigma_n (\omega)');
ylabel('\sigma_1(10^5\Omega^{-1}cm^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_Nb_real_cond_finite_tau_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');


% Plotting imaginary part -----------------------------------
figure
hold on
for i = 1:nT
    col = colors(i, :);

    % 1. BCS
    plot(x, sigma2(i, :)*sigma_dc, 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS $T = %.1f$~K', T(i)));
    % 2. BCS with gamma
    plot(x, sigma2_gamma(i, :)*sigma_dc, '--', 'Color', col, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('BCS-gamma $T = %.1f$~K', T(i)));
    % 3. Data
    idx_exp = abs(s2_data.T - T(i)) < 1e-6;
    if any(idx_exp)
        plot(s2_data.nu(idx_exp), s2_data.sigma(idx_exp), 'o', ...
            'DisplayName', sprintf('Data $T = %.1f$~K', T(i)), 'Color', col, 'MarkerFaceColor', col);
    end
end
grid on
title(sprintf('Imaginary part Conductivity of Nb ($T_c=%.1f~$K)', Tc), 'Interpreter', 'latex')
ylim([-0.1 8.5])
xlabel('Frequency (cm^{-1})'); 
%ylabel('\sigma_2/\sigma_n (\omega)');
ylabel('\sigma_2(10^5\Omega^{-1}cm^{-1})')
set(gca, 'FontSize', 18)
legend('Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
exportgraphics(gcf, 'BCS_Nb_imm_cond_finite_tau_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');