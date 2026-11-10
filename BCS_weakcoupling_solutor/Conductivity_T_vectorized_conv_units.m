
% Parameters
Tc     = 8.3;                           % K
kB     = 0.08617;                       % meV/K
alpha = 4.1;                           % coupling 2*D0/KbTc
Delta0 = 0.5*alpha*kB*Tc;               % meV, BCS weak coupling
c = physconst('LightSpeed') * 10^+2;    % c in cm
ratio_Tc = [0.001 0.3 0.6 0.8 0.99];
%T = ratio_Tc * Tc;
T = [4.5 6 7 7.5];                % Passing T in K directly
kT = T * kB;

f_THz = linspace(0.01, 1, 400);
hw = f_THz * 4.13567;
x = 10^12*f_THz/c;
nT = numel(T);
nW = numel(hw);

sigma1 = zeros(nT, nW);
sigma2 = zeros(nT, nW);

% Iterate over T
for i = 1:nT
    Delta = Delta0*tanh(1.74*sqrt(Tc/T(i) - 1));
    [s1, s2] = mattis_bardeen(hw, Delta, kT(i));
    sigma1(i, :) = s1; sigma2(i, :) = s2;
end

% Plotting real part --------------------------------------
figure
hold on
for i = 1:nT
    plot(x, sigma1(i, :), 'LineWidth', 1.5);
end
grid on
title(sprintf('Real part Conductivity of Nb ($T_c=%.1f~$K)',Tc ), 'Interpreter', 'latex')
ylim([-0.1 1.5])
xlabel('Frequency (cm^{-1})'); 
ylabel('\sigma_1/\sigma_n (\omega)');
set(gca, 'FontSize', 18)

%legend(arrayfun(@(r) sprintf('$T=%.1f\\,T_c$', r), ratio_Tc, 'UniformOutput', false), ...
    %'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')
legend(arrayfun(@(r) sprintf('$T=%.1f$~K', r), T, 'UniformOutput', false), ...
    'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')

exportgraphics(gcf, 'BCS_Nb_real_cond_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector');

% Plotting imaginary part -----------------------------------
figure
hold on
for i = 1:nT
    plot(x, sigma2(i, :), 'LineWidth', 1.5);
end
grid on
title(sprintf('Imaginary part Conductivity of Nb ($T_c=%.1f~$K)',Tc ), 'Interpreter', 'latex')
ylim([-0.1 4])
xlabel('Frequency (cm^{-1})'); 
ylabel('\sigma_2/\sigma_n (\omega)');
set(gca, 'FontSize', 18)

%legend(arrayfun(@(r) sprintf('$T=%.1f\\,T_c$', r), ratio_Tc, 'UniformOutput', false), ...
    %'Interpreter', 'latex', 'FontSize', 16, 'Location', 'best')
legend(arrayfun(@(r) sprintf('$T=%.1f$~K', r), T, 'UniformOutput', false), ...
   'Interpreter', 'latex', 'FontSize', 18, 'Location', 'best')

exportgraphics(gcf, 'BCS_Nb_imm_cond_diffT_conv_units_fitexp.pdf', 'ContentType', 'vector')