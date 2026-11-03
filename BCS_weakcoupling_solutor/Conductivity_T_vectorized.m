% Parameters
Tc     = 9.3;                 % K
kB     = 0.08617;             % meV/K
Delta0 = 0.5 * 3.65 * kB * Tc; % meV, BCS weak coupling

ratio_Tc = [0.01 0.3 0.6 0.8];
%T = ratio_Tc * Tc;
T = [4.5 6 7.5 8 8.2];
hw = linspace(0.01, 3.5, 400) * 2 * Delta0;
x = hw / (2 * Delta0);

nT = numel(T);
nW = numel(hw);

sigma1 = zeros(nT, nW);
sigma2 = zeros(nT, nW);

% Iterate over T
for i = 1:nT
    Delta = Delta0*tanh(1.74*sqrt(Tc/T(i) - 1));
    [s1, s2] = mattis_bardeen(hw, Delta, T(i));

    sigma1(i, :) = s1;
    sigma2(i, :) = s2;
end

% Plotting real part
figure
hold on
for i = 1:nT
    plot(x, sigma1(i, :), 'LineWidth', 1.5);
end
grid on
title('Real part Conductivity of Nb ($T_c=9.3~$K)', 'Interpreter', 'latex')
ylim([-0.1 2])
xlabel('$\hbar\omega / 2\Delta(0)$', 'Interpreter', 'latex');
ylabel('$\sigma/\sigma_n~(\omega)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 16)

%legend(arrayfun(@(r) sprintf('$T=%.1f\\,T_c$', r), ratio_Tc, 'UniformOutput', false), ...
    %'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')
legend(arrayfun(@(r) sprintf('$T=%.1f$~K', r), T, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')
print('BCS_Nb_real_cond_diffT_Maksimov.pdf', '-dpdf', '-bestfit')

% Plotting imaginary part
figure
hold on
for i = 1:nT
    plot(x, sigma2(i, :), 'LineWidth', 1.5);
end
grid on
title('Imaginary part Conductivity of Nb ($T_c=9.3~$K)', 'Interpreter', 'latex')
ylim([-0.1 2])
xlabel('$\hbar\omega / 2\Delta(0)$', 'Interpreter', 'latex');
ylabel('$\sigma/\sigma_n~(\omega)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 16)

legend(arrayfun(@(r) sprintf('$T=%.1f$~K', r), T, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')

print('BCS_Nb_imm_cond_diffT_Maksimov.pdf', '-dpdf', '-bestfit')
