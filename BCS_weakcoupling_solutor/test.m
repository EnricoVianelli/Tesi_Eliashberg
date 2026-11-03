% Parameters
Tc     = 9.3;               % K
kB     = 0.08617;           % meV/K
Delta0 = 0.5*3.65*kB*Tc;    % meV, BCS weak coupling
ratio_Tc = 0.8;
T = ratio_Tc*Tc;
Delta  = Delta0*tanh(1.74*sqrt(Tc/T - 1));    % approx. BCS gap(T) for T near Tc

% Setting up the integration
hw = linspace(0.01, 3.5, 400)*2*Delta0;        % Photon energies
[s1, s2] = mattis_bardeen(hw, Delta, T);       % Function to evaluate M-B integral conductivity

% Plotting
%title_text = ['Conductivity of Nb ($T_c=9.3~$K) at $T=$', num2str(ratio_Tc),'$T_c$'];
title_text = ['Conductivity of Nb ($T_c=9.3~$K) at $T=$', num2str(ratio_Tc),'$T_c$'];
file_name = ['BCS_Nb_conductivity_', num2str(ratio_Tc), 'Tc.pdf'];

figure
plot(hw/(2*Delta0), s1, 'r', hw/(2*Delta0), s2, 'g'); 
title(title_text, 'Interpreter', 'Latex')
ylim([-0.1 2])
xlabel('$\hbar\omega / 2\Delta(0)$', 'Interpreter', 'Latex'); 
ylabel('$\sigma/\sigma_n~(\omega)$', 'Interpreter', 'Latex');
legend('$\sigma_1/\sigma_n$','$\sigma_2/\sigma_n$', 'Interpreter', 'Latex', 'Fontsize', 16); grid on
set(gca, 'FontSize', 16) % Modifica il font di Assi, Titolo e Label
print(file_name, '-dpdf', '-bestfit')
