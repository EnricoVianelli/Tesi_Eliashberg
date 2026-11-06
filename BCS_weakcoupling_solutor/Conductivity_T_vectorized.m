% Parameters
Tc     = 7.2;                 % K
kB     = 0.08617;             % meV/K
alpha = 3.95;                 % coupling 2*D0/KbTc
Delta0 = 0.5*alpha*kB*Tc;     % meV, BCS weak coupling

ratio_Tc = [0.01 0.3 0.6 0.8];
T = ratio_Tc * Tc;
%T = [4.5 6 7.5 8 8.2];     % Passing T in K directly
kT = T * kB;
hw = linspace(0.01, 3.5, 400) * 2 * Delta0;
x = hw / (2 * Delta0);

nT = numel(T);
nW = numel(hw);

sigma1 = zeros(nT, nW);
sigma2 = zeros(nT, nW);

% Iterate over T
for i = 1:nT
    Delta = Delta0*tanh(1.74*sqrt(Tc/T(i) - 1));
    [s1, s2] = mattis_bardeen(hw, Delta, kT(i));

    sigma1(i, :) = s1;
    sigma2(i, :) = s2;
end

% Plotting real part --------------------------------------
figure
hold on
for i = 1:nT
    plot(x, sigma1(i, :), 'LineWidth', 1.5);
end
grid on
title(sprintf('Real part Conductivity of Pb ($T_c=%.1f~$K)',Tc ), 'Interpreter', 'latex')
ylim([-0.1 2])
xlabel('$\hbar\omega / 2\Delta(0)$', 'Interpreter', 'latex');
ylabel('$\sigma_1/\sigma_n~(\omega)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18)

legend(arrayfun(@(r) sprintf('$T=%.1f\\,T_c$', r), ratio_Tc, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
%legend(arrayfun(@(r) sprintf('$T=%.1f$~K', r), T, 'UniformOutput', false), ...
    %'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')

%print(, '-dpdf', '-bestfit')
exportgraphics(gcf, 'BCS_Pb_real_cond_diffT.pdf', 'ContentType', 'vector');

% Plotting imaginary part -----------------------------------
figure
hold on
for i = 1:nT
    plot(x, sigma2(i, :), 'LineWidth', 1.5);
end
grid on
title(sprintf('Imaginary part Conductivity of Pb ($T_c=%.1f~$K)',Tc ), 'Interpreter', 'latex')
ylim([-0.1 2])
xlabel('$\hbar\omega / 2\Delta(0)$', 'Interpreter', 'latex');
ylabel('$\sigma_2/\sigma_n~(\omega)$', 'Interpreter', 'latex');
set(gca, 'FontSize', 18)

legend(arrayfun(@(r) sprintf('$T=%.1f\\,T_c$', r), ratio_Tc, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')
%legend(arrayfun(@(r) sprintf('$T=%.1f$~K', r), T, 'UniformOutput', false), ...
%   'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')

%print('BCS_Pb_imm_cond_diffT.pdf', '-dpdf', '-bestfit')
exportgraphics(gcf, 'BCS_Pb_imm_cond_diffT.pdf', 'ContentType', 'vector')